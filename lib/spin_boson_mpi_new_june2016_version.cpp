/**
This code is to run the full version of the spin-boson simulation, using MPI.
Some of the data that was previously held in classes have been moved so that
they are easier to transfer around and manipulate.  Now, we will run parallel
simulations on each MPI unit (be it processor or node), with a GPU (or two?) 
working in tandem to calculate the action for parallel tempering moves. Each
unit will hold the lattice itself, its parameters (which will be exchanged in 
the parallel tempering), and the corresponding interactions.  The requirement
will be that the lattices are all the same size.  The main difference between
this and the version that uses the "LongRangeWolff" class is that it will allow
the user to define their own Monte Carlo step function, which makes it much more
flexible.
**/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "IsingLattice2D.h"
#include "class_mc_io.h"
#include "LongRangeWolff2D.h"
#include "Matrix.h"
#include "MemTimeTester.h"
#include "obs_calc_fast.cuh"
extern "C" {
#include "random.h"
}
#include <mpi.h>

using namespace std;

thrust::host_vector<double> IsingLattice2D::get_thrust_vector(){
	thrust::host_vector<double> result(Lx*Ly);
	for (int i = 0; i < Lx; ++i){
		for (int j = 0; j < Ly; ++j){
			result[i*Ly + j] = (double) spins[i][j];
		}
	}
	return result;
}

thrust::host_vector<double> GeneralLRW::get_thrust_interactions(){
	thrust::host_vector<double> result(interactions.size()*interactions[0].size());
	for (int i = 0; i < interactions.size(); ++i){
		for (int j = 0; j < interactions[i].size(); ++j){
			result[i*interactions[i].size() + j] = interactions[i][j];
		}
	}
	return result;
}

void ptemp(int num_procs, int id, double action, double alpha, double gamma, double sx, IsingLattice2D& lat){
    //use MPI to perform a parallel tempering step
    //when it comes time to send the lattices, master process will send first
    //and child processes will all receive first.  This might not be the quickest way
    //to do it but it ensures nothing gets locked up waiting to send/receive
    //
    //lattice message tags will be the receiving id

    //to begin, action = alpha*c + gamma*gamma_cont.  we want the variable "action" to be c, so action = (action - gamma*gamma_cont) / alpha
    MPI_Status Stat[2];
	MPI_Request req[2];
	std::vector<int> lat_buffer_out = lat.get_bool_spins();
	std::vector<int> lat_buffer_in(lat.get_N());
    double gamma_cont = ((double)lat.get_N())*(2.0*sx - 1);
    action = (action - gamma*gamma_cont) / alpha;

    //master process
    //1. receive action from other processes
    //2. test for switches, storing the send-to id in elements of an array
    //3. send the send-to id to each other process
    //4. send lattice if send-to[0] != 0
    if(id == 0){
        std::vector<double> actions(num_procs);
        std::vector<double> alphas(num_procs);
        std::vector<double> gamma_conts(num_procs);
        std::vector<int> receive_from(num_procs);
        std::vector<int> send_to(num_procs);
        actions[0] = action;
        alphas[0] = alpha;
        gamma_conts[0] = gamma_cont;
        for (int i = 1; i < num_procs; ++i){
            MPI_Recv(&(actions[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(alphas[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(gamma_conts[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            receive_from[i] = i;
        }

        //test for switches
        //assume for now that action = alpha*C + gamma*gamma_cont; this will need to be adjust later
        //start with lowest id and travel up
        double prob;
        for (int i = 0; i < num_procs - 1; ++i){
            prob = exp((alphas[receive_from[i]] - alphas[receive_from[i + 1]]) * (actions[receive_from[i]] - actions[receive_from[i + 1]]) + 2.0*gamma*(gamma_conts[receive_from[i]] - gamma_conts[receive_from[i + 1]]));
            //cout << "Parallel Tempering compare alphas " << alphas[receive_from[i]] << " and " << alphas[receive_from[i + 1]] << "; and actions " <<
            //            actions[receive_from[i]] << " and " << actions[receive_from[i + 1]] << "; and gamma contributions " <<
            //            gamma_conts[receive_from[i]] << " and " << gamma_conts[receive_from[i + 1]] << " for gamma = " << gamma << ", with probability " << prob << "\n";
            if (drand1_() < prob){
                receive_from[i + 1] = receive_from[i];
                receive_from[i] = i + 1;
            }
        }
        //invert receive_from to get send_to
        for (int i = 0; i < num_procs; ++i){
            send_to[receive_from[i]] = i;
        }
        //receive_from[i] now gives the id of the lattice that process i will receive
        //message each process and tell them which process they will send their lattice to and which process will receive their lattice
        //tags for send-to id will be the process id, tags for receive-from id will be num_procs + id
        for (int i = 1; i < num_procs; ++i){
            //send-to
            MPI_Send(&send_to[i], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            //receive-from
            MPI_Send(&receive_from[i], 1, MPI_DOUBLE, i, num_procs + i, MPI_COMM_WORLD);
        }
        //cout << "MPI master process: send_to = " << send_to[0] << ", receive_from = " << receive_from[0] << "\n";
	
        //send lat_buffer to send_to[0] and receive from receive_from[0]
        if(send_to[0] != 0){
            MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to[0], send_to[0], MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from[0], 0, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, Stat);
            lat.copy_bool_spins(lat_buffer_in);
        }
	}
    else{
        //child processes
        //1. send action to master process
        //2. receive send-to process and receive-from process
        //3. if send-to[id] != id, send lattice and receive new lattice
        MPI_Send(&action, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&alpha, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&gamma_cont, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);


        //cout << "MPI child process sending message action = " << action << ", alpha = " << alpha << "\n";

        int send_to, receive_from;
        MPI_Recv(&send_to, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &Stat[0]);
        MPI_Recv(&receive_from, 1, MPI_DOUBLE, 0, num_procs + id, MPI_COMM_WORLD, &Stat[0]);
        //cout << "MPI child process: send_to = " << send_to << ", receive_from = " << receive_from << "\n";

        //receive lattice and send lattice
        if(send_to != id){
            MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from, id, MPI_COMM_WORLD, &req[0]);
            //cout << "Process id " << id << " received lattice copy\n";
            MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to, send_to, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, Stat);
            lat.copy_bool_spins(lat_buffer_in);
        }
    }

    //cout << "Parallel tempering step completed\n";
}

int main(int argc, char* argv[]){
//
//  define global variables
//
    int p;//number of processes
    int id;//ID of this process
    int seed;//RNG seed for each process
    double random_value;//test random value
	bool gpu;//gpu, yay or nay, determined by command line argument

//
//  Initialize MPI.
//
    MPI_Init ( NULL, NULL );
//
//  Get the number of processes.
//
    MPI_Comm_size ( MPI_COMM_WORLD, &p );
//
//  Get the ID of this process.
//
    MPI_Comm_rank ( MPI_COMM_WORLD, &id );

	gpu = (argc <= 3);
	std::cout << "Process " << id << " GPU, yay or nay: " << gpu << "\n";

//
//  After MPI_Init is set, master node should do the following:
//  The master process prints a message.
//  Need a base parameter file: just specify the name of an input file with the parameters
//  Need to know which parameter to vary and how to vary it.  This should be compatible with # of processes
//
    class_mc_params params;
    ifstream infile;
    string vary_param;
    double delta_p, ind_var;
    infile.open(argv[1]);
    read_input_ising(&infile, &params);
    read_input_spin_boson(&infile, &(params.sbparams));
    read_input_mpi(&infile, &vary_param, &delta_p);
    //change parameter based on id, delta, and given param name.
    if (vary_param.compare("alpha") == 0){
        params.sbparams.A0 += id*delta_p;
        ind_var = params.sbparams.A0;
    }
    else if (vary_param.compare("delta") == 0){
        params.sbparams.delta += id*delta_p;
        ind_var = params.sbparams.delta;
    }
    else if (vary_param.compare("r") == 0){
        params.lengths[0] += id*delta_p;
        ind_var = params.lengths[0];
    }
    params.rand_seed += id*100;
    apply_spin_boson_params(&params);
    if (id == 0){
        cout << "Params for process " << id << ":\n" << params.to_string() << "\n\n\n";
    }

//
//  Simulation should be ready to set up for each process with different params and random seed
//  
    IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
    IsingLattice2D& latref = lat; //this will be the main workable item, staying local to its original process in the node
    GeneralLRW wolff = GeneralLRW(params);
    wolff.set_mag(latref);
    //cout << "Interactions and Interaction sum for process " << id << ":\n\n";
    //wolff.print_interactions();
    //wolff.print_interaction_sum();
    //wolff.print_site_sums();

//
//  Set up measurement apparatus
//
    class_mc_measurements results;
    results.names = {"loc", "loc2", "loc4", "loc_bind", "mag", "mag2", "mag4", "mag_bind",  "x_corr", "sx", "cluster"};
    results.values = {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}};
    MemTimeTester timer;
    timer.flag_start_time("full simulation");
    double loc, mag;


    if (gpu){
        //
        //  Define parallel tempering thrust variables
        //
        thrust::host_vector<double> thrustlat = latref.get_thrust_vector();
        thrust::host_vector<double>& thrustlat_ref = thrustlat;
        thrust::host_vector<double> thrustint = wolff.get_thrust_interactions();
        thrust::host_vector<double>& thrustint_ref = thrustint;
        double fast_action;

        
        //
        //  Equilibrate, using parallel tempering
        //
        timer.flag_start_time("equilibration");
        for (int i = 0; i < (params.max_dumps/10); ++i){
            for (int j = 0; j < params.eq_time/(params.max_dumps/10); ++j){
                wolff.step(latref);
            }
            thrustlat = latref.get_thrust_vector();
            timer.flag_start_time("action calculation");
            fast_action = cufft_calc_action(thrustlat_ref, thrustint_ref, latref.get_Lx(), latref.get_Ly());
            timer.flag_end_time("action calculation");
            timer.flag_start_time("parallel tempering");
            ptemp(p, id, fast_action, params.sbparams.A0, params.Js[1], wolff.calc_sx(latref), latref);
            wolff.set_mag(latref);
            timer.flag_end_time("parallel tempering");
        }


        //
        //  Run the measurement portion of the simulation
        //  Use each "dump" as the parallel tempering step
        //
        //dim3 threads(lat.get_Ly(), 1, 1);
        for (int dump = 0; dump < params.max_dumps; ++dump){
            for (int measure = 0; measure < params.measures_per_dump; ++measure){
                timer.flag_start_time("steps");
                    for (int n = 0; n < params.steps_per_measure; ++n){
                        wolff.step(latref);
                    }
                    timer.flag_start_time("measurements");
                    loc = wolff.calc_loc(latref);
                    mag = wolff.get_mag();

                    results.record("loc", loc);
                    results.record("loc2", loc*loc);
                    results.record("loc4", loc*loc*loc*loc);

                    results.record("mag", mag);
                    results.record("mag2", mag*mag);
                    results.record("mag4", mag*mag*mag*mag);

                    results.record("cluster", ((double)wolff.get_cluster_size()) / lat.get_N());
                    results.record("x_corr", wolff.calc_s1s2(latref));
                    results.record("sx", wolff.calc_sx(latref));
                    timer.flag_end_time("measurements");

                timer.flag_end_time("steps");
            }
            thrustlat = latref.get_thrust_vector();
            //parallel tempering
            timer.flag_start_time("action calculation");
            //fast_action = thrust_calc_action_general(thrustlat_ref, thrustint_ref, latref.get_Lx(), latref.get_Ly(), threads);
            fast_action = cufft_calc_action(thrustlat_ref, thrustint_ref, latref.get_Lx(), latref.get_Ly());
            timer.flag_end_time("action calculation");
            timer.flag_start_time("parallel tempering");
            ptemp(p, id, fast_action, params.sbparams.A0, params.Js[1], wolff.calc_sx(latref), latref);
            wolff.set_mag(latref);
            timer.flag_end_time("parallel tempering");
        }
    }
    else{

        //
        //  Equilibrate
        //
        timer.flag_start_time("equilibration");
        for (int i = 0; i < params.eq_time; ++i){
            wolff.step(latref);
        }
        timer.flag_end_time("equilibration");

        for (int dump = 0; dump < params.max_dumps; ++dump){
            for (int measure = 0; measure < params.measures_per_dump; ++measure){
                timer.flag_start_time("steps");
                for (int n = 0; n < params.steps_per_measure; ++n){
                    wolff.step(latref);
                }
                timer.flag_start_time("measurements");
                loc = wolff.calc_loc(latref);
                mag = wolff.get_mag();

                results.record("loc", loc);
                results.record("loc2", loc*loc);
                results.record("loc4", loc*loc*loc*loc);

                results.record("mag", mag);
                results.record("mag2", mag*mag);
                results.record("mag4", mag*mag*mag*mag);

                results.record("cluster", ((double)wolff.get_cluster_size()) / lat.get_N());
                results.record("x_corr", wolff.calc_s1s2(latref));
                results.record("sx", wolff.calc_sx(latref));

                timer.flag_end_time("measurements");
                timer.flag_end_time("steps");
            }
        }
    }
//
//  test fast action calculation
//
/*
    double slow_action = wolff.calc_action_slow(latref);
    fast_action = thrust_calc_action_general(thrustlat_ref, thrustint_ref, latref.get_Lx(), latref.get_Ly());
    std::stringstream outstring;
    outstring << "Printing action test for process " << id << "\nLattice:\n" << latref.to_string() << "\n" << wolff.get_int_string() << "Slow action: " << slow_action << ", Fast action: " << fast_action << "\n";
    cout << outstring.str();
*/


//
//  Write the results
//
    double loc_avg, loc2_avg, loc4_avg, loc_bind_avg, loc_bind_err, mag_avg, mag2_avg, mag4_avg, mag_bind_err, mag_bind_avg, xcorr_avg, sx_avg, cluster_avg;
    loc_avg = mean(results.get_vals("loc"));
    loc2_avg = mean(results.get_vals("loc2"));
    loc4_avg = mean(results.get_vals("loc4"));
    loc_bind_avg = loc4_avg/loc2_avg/loc2_avg;
    loc_bind_err = bootstrap(results.get_vals("loc2"));

    mag_avg = mean(results.get_vals("mag"));
    mag2_avg = mean(results.get_vals("mag2"));
    mag4_avg = mean(results.get_vals("mag4"));
    mag_bind_avg = mag4_avg/mag2_avg/mag2_avg;
    mag_bind_err = bootstrap(results.get_vals("mag2"));

    xcorr_avg = mean(results.get_vals("x_corr"));
    sx_avg = mean(results.get_vals("sx"));
    cluster_avg = mean(results.get_vals("cluster"));

    cout << "Results for process #" << id << ":\n";
    timer.flag_end_time("full simulation");
    timer.print_timers();
    cout << "total steps = " << params.steps_per_measure*params.measures_per_dump*params.max_dumps << ", parallel tempering steps = " << params.max_dumps << "\n";
    write_outputs_var(id, results);//write measurements to a dump file

//
//  Send results to master process for easy writing to a file
//
    if(id == 0){
        std::vector<double> ind_vars(p);//vector of independent variables (can change from alpha if this changes in earlier parts of the code)
        std::vector<double> locs(p), loc2s(p), loc_binds(p), loc_bind_errs(p), mags(p), mag2s(p), mag_binds(p), mag_bind_errs(p), xcorrs(p), sxs(p), clusters(p);
        ind_vars[0] = ind_var;

        locs[0] = loc_avg;
        loc2s[0] = loc2_avg;
        loc_binds[0] = loc_bind_avg;
        loc_bind_errs[0] = loc_bind_err;

        mags[0] = mag_avg;
        mag2s[0] = mag2_avg;
        mag_binds[0] = mag_bind_avg;
        mag_bind_errs[0] = mag_bind_err;

        xcorrs[0] = xcorr_avg;
        sxs[0] = sx_avg;
        clusters[0] = cluster_avg;

        MPI_Status Stat;
        for (int i = 1; i < p; ++i){
            MPI_Recv(&(ind_vars[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(locs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc2s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc_binds[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc_bind_errs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(mags[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag2s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag_binds[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag_bind_errs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(xcorrs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(sxs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(clusters[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

        }
        //write results
        FILE * file;
        file = fopen("results.dat", "w");
	    fprintf(file, "%-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s\n", 
                vary_param, "loc", "loc2", "loc_bind", "lb_err", "mag", "mag2", "mag_bind", "mb_err", "x_corr", "sx", "cluster");
        for (int i = 0; i < p; ++i){
            fprintf(file, "%-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n", 
                    ind_var[i], locs[i], loc2s[i], loc_binds[i], loc_bind_errs[i], mags[i], mag2s[i], mag_binds[i], mag_bind_errs[i], xcorrs[i], sxs[i], clusters[i]);
        }
        fclose (file);
    }
    else{
        MPI_Send(&ind_var, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&loc_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc2_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc_bind_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc_bind_err, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&mag_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag2_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag_bind_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag_bind_err, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&xcorr_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&sx_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&cluster_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        
    }
	
//
//  End the MPI process
//
    MPI_Finalize();

return 0;
}
