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
#include <thrust/inner_product.h>
extern "C" {
#include "random.h"
}
#include <mpi.h>
//#include <helper_cuda.h>

using namespace std;

double ptemp(int num_procs, int id, double action, double alpha, double gamma, double sx, double J, double sflips, IsingLattice2D& lat, bool *switched){
    //use MPI to perform a parallel tempering step
    //when it comes time to send the lattices, master process will send first
    //and child processes will all receive first.  This might not be the quickest way
    //to do it but it ensures nothing gets locked up waiting to send/receive
    //
    //lattice message tags will be the receiving id

    //return value is the new action value

    //Don't use this function if alpha=0

    //"switched" tells if the lattice is moved or not
    *switched = false;//lattice does not move

    //to begin, action = alpha*c + gamma*gamma_cont + J*J_cont.  we want the variable "action" to be c, so action = (action - gamma*gamma_cont - J*J_cont) / alpha
    MPI_Status Stat[2];
	MPI_Request req[2];
	std::vector<int> lat_buffer_out = lat.get_bool_spins();
    lat_buffer_out.shrink_to_fit();
	std::vector<int> lat_buffer_in(lat.get_N());
    lat_buffer_in.shrink_to_fit();
    double gamma_cont = ((double)lat.get_N())*(2.0*sx - 1);
    double J_cont = ((double)lat.get_N())*(2.0*sflips - 1);
    double new_action;
    action = (action - gamma*gamma_cont - J*J_cont) / alpha;

    //master process
    //1. receive action from other processes
    //2. test for switches, storing the send-to id in elements of an array
    //3. send the send-to id to each other process
    //4. send lattice if send-to[0] != 0
    if(id == 0){
        std::vector<double> actions(num_procs);
        std::vector<double> alphas(num_procs);
        std::vector<double> gamma_conts(num_procs);
        std::vector<double> J_conts(num_procs);
        std::vector<double> new_action_list(num_procs);
        std::vector<int> receive_from(num_procs);
        std::vector<int> send_to(num_procs);
        actions[0] = action;
        alphas[0] = alpha;
        gamma_conts[0] = gamma_cont;
        J_conts[0] = J_cont;
        receive_from[0] = 0;
        send_to[0] = 0;
        for (int i = 1; i < num_procs; ++i){
            MPI_Recv(&(actions[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(alphas[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(gamma_conts[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(J_conts[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            receive_from[i] = i;
        }

        //test for switches
        //assume for now that action = alpha*C + gamma*gamma_cont
        //start with lowest id and travel up
        for (int i = 0; i < num_procs - 1; ++i){
            double prob;
            prob = exp((alphas[receive_from[i]] - alphas[receive_from[i + 1]]) * (actions[receive_from[i]] - actions[receive_from[i + 1]]));
            //std::cout << "Action difference: " << actions[receive_from[i]] - actions[receive_from[i + 1]] << "\n";
            //         + 2.0*gamma*(gamma_conts[receive_from[i]] - gamma_conts[receive_from[i + 1]]));this part is probably wrong

            if (drand1_() < prob){
                //if prob > 1, that means the switching i and i + 1 has a favorable change in action
                //if prob < 1, the switch has an unfavorable change in action but there is still some probability of switching
                receive_from[i + 1] = receive_from[i];
                receive_from[i] = i + 1;
            }
        }
        //invert receive_from to get send_to
        //list new actions for the lattices each process will be receiving
        for (int i = 0; i < num_procs; ++i){
            send_to[receive_from[i]] = i;
            new_action_list[i] = alphas[i] * actions[receive_from[i]] + gamma * gamma_conts[receive_from[i]] + J * J_conts[receive_from[i]];
        }
        new_action = new_action_list[0];
        //receive_from[i] now gives the id of the lattice that process i will receive
        //message each process and tell them which process they will send their lattice to and which process will receive their lattice
        //tags for send-to id will be the process id, tags for receive-from id will be num_procs + id, tags for new actions will be num_procs + 2*id
        for (int i = 1; i < num_procs; ++i){
            //send-to
            MPI_Send(&(send_to[i]), 1, MPI_INT, i, i, MPI_COMM_WORLD);
            //receive-from
            MPI_Send(&(receive_from[i]), 1, MPI_INT, i, num_procs + i, MPI_COMM_WORLD);
            //send new action
            MPI_Send(&(new_action_list[i]), 1, MPI_DOUBLE, i, num_procs + 2*i, MPI_COMM_WORLD);
        }
        //cout << "MPI master process: send_to = " << send_to[0] << ", receive_from = " << receive_from[0] << "\n";
	
        //send lat_buffer to send_to[0] and receive from receive_from[0]
        if(send_to[0] != 0){
            MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to[0], send_to[0], MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from[0], 0, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, Stat);
            lat.copy_bool_spins(lat_buffer_in);
            *switched = true; //lattice is moved
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
        MPI_Send(&J_cont, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);


        //cout << "MPI child process sending message action = " << action << ", alpha = " << alpha << "\n";

        int send_to, receive_from;
        MPI_Recv(&send_to, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &Stat[0]);
        MPI_Recv(&receive_from, 1, MPI_INT, 0, num_procs + id, MPI_COMM_WORLD, &Stat[0]);
        MPI_Recv(&new_action, 1, MPI_DOUBLE, 0, num_procs + 2*id, MPI_COMM_WORLD, &Stat[0]);
        //cout << "MPI child process: send_to = " << send_to << ", receive_from = " << receive_from << "\n";

        //receive lattice and send lattice
        if(send_to != id){
            MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from, id, MPI_COMM_WORLD, &req[0]);
            //cout << "Process id " << id << " received lattice copy\n";
            MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to, send_to, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, Stat);
            lat.copy_bool_spins(lat_buffer_in);
            *switched = true; //lattice is moved
        }

    }

    return new_action;
}

int main(int argc, char* argv[]){
//
//  define global variables
//
    int p;//number of processes
    int id;//ID of this process
    int seed;//RNG seed for each process
    int device_count = 0;//number of GPUs
    int proc_per_device;//number of threads per GPU
    double random_value;//test random value
	bool gpu;//gpu, yay or nay, determined by command line argument
    bool mpi_use;//mpi, yay or nay, determined by command line argument
    cudaDeviceProp gpu_stats;
    std::stringstream outstring;

    mpi_use = (argc <= 3);
    if(mpi_use){
//
//  Initialize MPI.
//
    MPI_Init ( &argc, &argv );
//
//  Get the number of processes.
//
    MPI_Comm_size ( MPI_COMM_WORLD, &p );
//
//  Get the ID of this process.
//
    MPI_Comm_rank ( MPI_COMM_WORLD, &id );
    }
    else{
        p = 1;
        id = 0;
        std::cout << "Not using MPI\n";
    }

	gpu = (argc <= 4);
	if(id == 0){
        std::cout << "Using GPU, yes or no: " << (gpu ? "yes" : "no") << "\n";
    }


//
//  Check GPU properties
//
    if(gpu){
        cudaGetDeviceCount(&device_count);
        proc_per_device = p / device_count;
        if(id/proc_per_device < device_count){
            cudaSetDevice(id/proc_per_device);
        }
        else{
            cudaSetDevice(0);
        }
        if(id == 0){
            for(int i = 0; i < p; ++i){
                outstring << "Thread " << i << " set to Device " << i/proc_per_device << "\n";
            }
        }
        for(int d = 0; d < device_count; ++d){
            if(id == d * proc_per_device){
                cudaGetDeviceProperties(&gpu_stats, d);
                outstring << "GPU properties:\nDevice: " << d << "\nName: " << gpu_stats.name << "\nThreads per block: " << gpu_stats.maxThreadsPerBlock
                    << "\nThread dimensions: {" << gpu_stats.maxThreadsDim[0] << "," << gpu_stats.maxThreadsDim[1] << "," << gpu_stats.maxThreadsDim[2]
                    << "}\nMax grid size: {" << gpu_stats.maxGridSize[0] << "," << gpu_stats.maxGridSize[1] << "," << gpu_stats.maxGridSize[2] << "}\n\n";
                std::cout << outstring.str();
                outstring.str("");
                show_memory();
            }
        }
    }


//
//  After MPI_Init is set, master node should do the following:
//  The master process prints a message.
//  Need a base parameter file: just specify the name of an input file with the parameters
//  Need to know which parameter to vary and how to vary it.  This should be compatible with # of processes
//
    class_mc_params params;
    ifstream infile;
    string vary_param;
    double ind_var;
    std::vector<double> param_vals;
    std::vector<double>& param_ref = param_vals;
    infile.open(argv[1]);
    read_input_ising(&infile, &params);
    read_input_spin_boson(&infile, &(params.sbparams));
    read_input_mpi(&infile, &vary_param, param_ref);
    if(param_ref.size() != p && id == 0){
        std::cout << "Error reading parallel tempering params! Number of param values in input file does not match number of processes\n";
    }
    //change parameter based on id, delta, and given param name.
    if (vary_param.compare("alpha") == 0){
        params.sbparams.A0 = param_ref[id];
        ind_var = params.sbparams.A0;
    }
    else if (vary_param.compare("J") == 0){
        params.Js[0] = param_ref[id];
        ind_var = params.Js[0];
    }
    else if (vary_param.compare("delta") == 0){
        params.sbparams.delta = param_ref[id];
        ind_var = params.sbparams.delta;
    }
    else if (vary_param.compare("r") == 0){
        params.spacings[0] = param_ref[id];
        ind_var = params.spacings[0];
    }
    else if (vary_param.compare("Nx") == 0){
        params.lengths[0] = param_ref[id];
        ind_var = params.lengths[0];
    }
    else if (vary_param.compare("Nt") == 0){
        params.lengths[1] = param_ref[id];
        ind_var = params.lengths[1];
    }
    else if (vary_param.compare("beta") == 0){
        params.beta = param_ref[id];
        ind_var = params.beta;
    }
    params.rand_seed += id*100;
    apply_spin_boson_params(&params);
    if (id == 0){
        outstring << "Params for process " << id << ":\n" << params.to_string() << "\n\n\n";
        cout << outstring.str();
        makePath("./dump");
        outstring.str("");
    }


//
//  Simulation should be ready to set up for each process with different params and random seed
//  
    MemTimeTester timer;
    timer.flag_start_time("full simulation");
    IsingLattice2D lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
    IsingLattice2D& latref = lat; //this will be the main workable item, staying local to its original process in the node
    if(params.lattice.compare("ferro") == 0){
        lat.set_ferro(1);
    }
    else if(params.lattice.compare("antiferro") == 0){
        lat.set_antiferro();
    }
    else if(params.lattice.compare("stripe_x") == 0){
        lat.set_stripe_x();
    }
    else if(params.lattice.compare("stripe_y") == 0){
        lat.set_stripe_y();
    }
    else if(params.lattice.compare("given") == 0){
        read_state_pbm(id, latref);
    }
    GeneralLRW wolff = GeneralLRW(params);
    if(params.lattice.compare("vals") == 0){
        lat.set_signs(wolff.get_interactions());
    }
    //std::cout << "successfully created wolff object\n";
    wolff.set_mag(latref);    

    //Write interactions
    char intfile_name[100];
    sprintf(intfile_name, "dump/interactions%d.csv", id);
    std::ofstream file;
    file.open(intfile_name);
    file << lat.get_Lx() << "," << lat.get_Ly() << "\n";
    file << vec2str(wolff.get_interaction_vector()) << "\n";
    file.close();

//
//  Set up measurement apparatus
//
    class_mc_measurements results;
    results.names = {"loc", "loc2", "loc4", "xmag", "xmag2", "xmag4", "mag", "mag2", "mag4", "sx", "xkinks", "stagger_mag", "action", "cluster"};
    results.values = {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}};
    results.keep_function_names = {"space_corr"};
    results.keep_function_vals = {{{}}};
    results.func_names = {"corr", "vort_corr"};
    results.function_num_measures = {0, 0};
    results.functions = {{}, {}};
    double loc, xmag, mag, ptemp_moves = 0, ptemp_total = 0, traversal = 0;
    int num_meas = 0, num_accept = 0, num_step = 0;

    if (gpu && params.sbparams.A0 != 0){
        //
        //  Define thrust variables
        //

        //for passing into LRW step
        thrust::host_vector<double> thrustsite_trash(latref.get_Ly()), thrustint_site_trash(latref.get_Ly()), site_corr_trash(latref.get_Ly());
        thrust::host_vector<double> thrustlat_trash(latref.get_Lx()*latref.get_Ly()), thrustint_trash(latref.get_Lx()*latref.get_Ly()), full_corr_trash(latref.get_Lx()*latref.get_Ly());
        thrust::host_vector<double>& thrustsite_trash_ref = thrustsite_trash;
        thrust::host_vector<double>& thrustint_site_trash_ref = thrustint_site_trash;
        thrust::host_vector<double>& site_corr_trash_ref = site_corr_trash;
        thrust::host_vector<double>& thrustlat_trash_ref = thrustlat_trash;
		thrust::host_vector<double>& thrustint_trash_ref = thrustint_trash;
        thrust::host_vector<double>& full_corr_trash_ref = full_corr_trash;

        //for data and calculations
        thrust::host_vector<double> thrustlat;
        thrustlat.shrink_to_fit();
        thrust::host_vector<double>& thrustlat_ref = thrustlat;
        lat.get_thrust_vector(thrustlat_ref);

        thrust::host_vector<double> thrustint;
        thrustint.shrink_to_fit();
        thrust::host_vector<double>& thrustint_ref = thrustint;
        wolff.get_thrust_interactions(thrustint_ref);

        thrust::host_vector<double> correlation(latref.get_Lx()*latref.get_Ly());
        correlation.shrink_to_fit();
        thrust::host_vector<double>& temp_corr = correlation;
        std::vector<double> corr_measure(latref.get_Lx()*latref.get_Ly());
        corr_measure.shrink_to_fit();

        thrust::host_vector<double> vort_state(latref.get_Lx()*latref.get_Ly());
        vort_state.shrink_to_fit();
        thrust::host_vector<double> vort_corr(latref.get_Lx()*latref.get_Ly());
        vort_corr.shrink_to_fit();
        thrust::host_vector<double>& vort_state_ref = vort_state;
        thrust::host_vector<double>& vort_corr_ref = vort_corr;
        std::vector<double> vort_corr_measure(latref.get_Lx()*latref.get_Ly());
        vort_corr_measure.shrink_to_fit();

        std::vector<double> space_corr_measure(latref.get_Lx());
        space_corr_measure.shrink_to_fit();
        
        double fast_action = 0, metropolis_action = 0;
        bool ptemp_switch = false;

        //
        //  Initialize CUFFT workspace
        //

        cufftHandle full_forward_plan, full_backward_plan, onesite_forward_plan, onesite_backward_plan;
        cudaError cuda_status;
        cufftDoubleReal *full_state_rs, *onesite_state_rs;//real space state
        cufftDoubleComplex *full_state_ft, *onesite_state_ft;//fourier space state
        int Lx = lat.get_Lx(), Ly = lat.get_Ly();
        int n[2] = {Lx, Ly};
        int n1[2] = {1, Ly};
        cuda_status = cudaMalloc((void**)&full_state_rs, sizeof(cufftDoubleReal)*Lx*Ly);
        if(cuda_status != cudaSuccess){
            outstring << "Process " << id << " cudaMalloc failed, " << cudaGetErrorString(cuda_status) << "; On device " << id/proc_per_device << "\n";
            std::cout << outstring.str();
            outstring.str("");
            show_memory();
        }
        cuda_status = cudaMalloc((void**)&full_state_ft, sizeof(cufftDoubleComplex)*Lx*(Ly / 2 + 1));
        if(cuda_status != cudaSuccess){
            outstring << "Process " << id << " cudaMalloc failed, " << cudaGetErrorString(cuda_status) << "; On device " << id/proc_per_device << "\n";
            std::cout << outstring.str();
            outstring.str("");
            show_memory();
        }
        cuda_status = cudaMalloc((void**)&onesite_state_rs, sizeof(cufftDoubleReal)*Ly);
        if(cuda_status != cudaSuccess){
            outstring << "Process " << id << " cudaMalloc failed, " << cudaGetErrorString(cuda_status) << "; On device " << id/proc_per_device << "\n";
            std::cout << outstring.str();
            outstring.str("");
            show_memory();
        }
        cuda_status = cudaMalloc((void**)&onesite_state_ft, sizeof(cufftDoubleComplex)*(Ly / 2 + 1));
        if(cuda_status != cudaSuccess){
            outstring << "Process " << id << " cudaMalloc failed, " << cudaGetErrorString(cuda_status) << "; On device " << id/proc_per_device << "\n";
            std::cout << outstring.str();
            outstring.str("");
            show_memory();
        }
        if(id == 0){
            cout << "Creating full forward plan\n";
        }
        if(cufftPlanMany(&full_forward_plan, 2, n, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, 1) != CUFFT_SUCCESS) {
            fprintf(stderr, "CUFFT Error: Unable to create plan\n");
            return 0;
        }
        if(id == 0){
            cout << "Creating full backward plan\n";
        }
        if (cufftPlanMany(&full_backward_plan, 2, n, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, 1) != CUFFT_SUCCESS) {
            fprintf(stderr, "CUFFT Error: Unable to create plan\n");
            return 0;
        }
        if(id == 0){
            cout << "Creating onesite forward plan\n";
        }
        if(cufftPlanMany(&onesite_forward_plan, 2, n1, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, 1) != CUFFT_SUCCESS) {
            fprintf(stderr, "CUFFT Error: Unable to create plan\n");
            return 0;
        }
        if(id == 0){
            cout << "Creating onesite backward plan\n";
        }
        if (cufftPlanMany(&onesite_backward_plan, 2, n1, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, 1) != CUFFT_SUCCESS) {
            fprintf(stderr, "CUFFT Error: Unable to create plan\n");
            return 0;
        }

        if(id == 0){ std::cout << "Starting Simulation\n\n";}
        lat.get_thrust_vector(thrustlat_ref);
        fast_action = cufft_calc_action_prealloc(thrustlat_ref, thrustint_ref, temp_corr, lat.get_Lx(), lat.get_Ly(), 
						full_forward_plan, full_backward_plan, cuda_status, full_state_rs, full_state_ft);
        metropolis_action = fast_action;

        //  
        //  Run Measurements
        //  At each measurement iteration, run the steps, take or throw away the measurement, then do parallel tempering
        //
        for(int measure = 1; measure <= (params.kept_measures + params.throwaway_measures); ++measure){
            
            timer.flag_start_time("step calculation");
            //Metropolis
            for (int n = 0; n < params.metropolis_steps; ++n){
                if(wolff.metropolis_step(latref, &metropolis_action)){
                    ++num_accept;
                    traversal += 1.0/params.lengths[0]/params.lengths[1];
                }
                ++num_step;
            }

            //Onesite
            for (int n = 0; n < params.onesite_steps; ++n){
                if(wolff.step_one_site_prealloc_lowmem(latref, &fast_action, temp_corr, thrustsite_trash_ref,
                    thrustint_site_trash_ref,site_corr_trash_ref,thrustlat_trash_ref,thrustint_trash_ref,full_corr_trash_ref,
                    &full_forward_plan, &full_backward_plan, &onesite_forward_plan, &onesite_backward_plan,
                    cuda_status, full_state_rs, full_state_ft, onesite_state_rs, onesite_state_ft)){
                        ++num_accept;
                        traversal += min(wolff.get_cluster_size()/params.lengths[0]/params.lengths[1], 1 - wolff.get_cluster_size()/params.lengths[0]/params.lengths[1]);
                }
                ++num_step;
            }

            //Wolff
            for (int n = 0; n < params.wolff_steps; ++n){
                wolff.step(latref);
                traversal += min(wolff.get_cluster_size()/params.lengths[0]/params.lengths[1], 1 - wolff.get_cluster_size()/params.lengths[0]/params.lengths[1]);
            }
            timer.flag_end_time("step calculation");

            //Measure
            timer.flag_start_time("measurements");
            lat.get_thrust_vector(thrustlat_ref);
            fast_action = cufft_calc_action_prealloc(thrustlat_ref, thrustint_ref, temp_corr, lat.get_Lx(), lat.get_Ly(), 
                full_forward_plan, full_backward_plan, cuda_status, full_state_rs, full_state_ft);
            wolff.set_mag(latref);
            if(measure > params.throwaway_measures){
                loc = wolff.calc_loc(latref);
                xmag = wolff.calc_xmag(latref);
                mag = wolff.calc_mag(latref);

                results.record("loc", loc);
                results.record("loc2", loc*loc);
                results.record("loc4", loc*loc*loc*loc);

                results.record("xmag", xmag);
                results.record("xmag2", xmag*xmag);
                results.record("xmag4", xmag*xmag*xmag*xmag);

                results.record("mag", mag);
                results.record("mag2", mag*mag);
                results.record("mag4", mag*mag*mag*mag);

                results.record("sx", wolff.calc_sx(latref));
                results.record("xkinks", wolff.calc_space_kinks(latref));
                results.record("stagger_mag", wolff.calc_sz_stagger(latref));
                results.record("action", fast_action);
                results.record("cluster", wolff.get_cluster_size());
                
                lat.get_vort_state(vort_state_ref);
                cufft_calc_action_prealloc(vort_state_ref, thrustint_ref, vort_corr_ref, lat.get_Lx(), lat.get_Ly(), 
                    full_forward_plan, full_backward_plan, cuda_status, full_state_rs, full_state_ft); //calculate vortex correlation
                for(int i = 0; i < temp_corr.size(); ++i){
                    corr_measure[i] = temp_corr[i];
                    vort_corr_measure[i] = vort_corr_ref[i];
                }
                for(int x = 0; x < lat.get_Lx(); ++x){
                    space_corr_measure[x] = 0.0;
                    for(int y = 0; y < lat.get_Ly(); ++y){
                        space_corr_measure[x] += corr_measure[x*lat.get_Ly() + y];
                    }
                }
                
                results.record_keep_function("space_corr", space_corr_measure);
                results.record("corr", corr_measure);
                results.record("vort_corr", vort_corr_measure);
                ++num_meas;
            }
            timer.flag_end_time("measurements");

            //parallel tempering
            if(mpi_use && measure%params.ptemp_steps == 0){
                //shuffle the lattices
                timer.flag_start_time("parallel tempering");
                fast_action = ptemp(p, id, fast_action, params.sbparams.A0, params.Js[1], wolff.calc_sx(latref), 
                                    params.Js[0], wolff.calc_space_kinks(latref), latref, &ptemp_switch);
                if (ptemp_switch){
                    ptemp_moves += 1.0;
                    traversal += 1.0;
                }
                wolff.set_mag(latref);
                timer.flag_end_time("parallel tempering");
                ptemp_total += 1.0;
            }

            
            if(id == 0){
                std::cout << "measurement step " << measure << " out of " 
                    << (params.kept_measures + params.throwaway_measures) << " completed\n";
            }
        }

        lat.get_thrust_vector(thrustlat_ref);
        if(abs(fast_action - cufft_calc_action(thrustlat_ref, thrustint_ref, temp_corr, latref.get_Lx(), latref.get_Ly())) >= 0.01){
            std::stringstream message;
            message << "Warning: fast action and cufft calculation do not match!  fast_action: " 
                    << fast_action << "; current cufft action: " << cufft_calc_action(thrustlat_ref, thrustint_ref, temp_corr, latref.get_Lx(), latref.get_Ly()) << "\n";
            std::cout << message.str();
        }

        //test to see if gpu memory availability has been depleted
        for(int d = 0; d < device_count; ++d){
            if(id == d * proc_per_device){
                outstring << "Memory usage for device " << d << ":\n";
                std::cout << outstring.str();
                outstring.str("");
                show_memory();
            }
        }
    }
    else{

        if(id == 0){
            std::cout << "Can only use wolff steps and no parallel tempering\n";
        }
        std::vector<double> corr_dummy(1);

        //  
        //  Run Measurements
        //  At each measurement iteration, run the steps, take or throw away the measurement, then do parallel tempering
        //
        for(int measure = 1; measure <= (params.kept_measures + params.throwaway_measures); ++measure){
            
            timer.flag_start_time("step calculation");

            //Wolff
            for (int n = 0; n < params.wolff_steps; ++n){
                wolff.step(latref);
                traversal += min(wolff.get_cluster_size()/params.lengths[0]/params.lengths[1], 1 - wolff.get_cluster_size()/params.lengths[0]/params.lengths[1]);
            }
            timer.flag_end_time("step calculation");

            //Measure
            timer.flag_start_time("measurements");
            wolff.set_mag(latref);
            if(measure > params.throwaway_measures){
                loc = wolff.calc_loc(latref);
                xmag = wolff.calc_xmag(latref);
                mag = wolff.calc_mag(latref);

                results.record("loc", loc);
                results.record("loc2", loc*loc);
                results.record("loc4", loc*loc*loc*loc);

                results.record("xmag", xmag);
                results.record("xmag2", xmag*xmag);
                results.record("xmag4", xmag*xmag*xmag*xmag);

                results.record("mag", mag);
                results.record("mag2", mag*mag);
                results.record("mag4", mag*mag*mag*mag);

                results.record("sx", wolff.calc_sx(latref));
                results.record("xkinks", wolff.calc_space_kinks(latref));
                results.record("stagger_mag", wolff.calc_sz_stagger(latref));
                results.record("action", 1.0);
                results.record("cluster", wolff.get_cluster_size());
                
                
                results.record("corr", corr_dummy);
                results.record("vort_corr", corr_dummy);
                ++num_meas;
            }
            timer.flag_end_time("measurements");
            std::cout << "measurement step " << measure << " out of " 
                << (params.kept_measures + params.throwaway_measures) << " completed\n";
        }
    }

//
//  Write the results
//
    double loc_avg, loc2_avg, loc4_avg, loc_bind_avg, loc_bind_err, 
            xmag_avg, xmag2_avg, xmag4_avg, xmag_bind_err, xmag_bind_avg, 
            mag_avg, mag2_avg, mag4_avg, mag_bind_err, mag_bind_avg, 
            action_avg, sx_avg, xkinks_avg, stagger_mag_avg, cluster_avg, ptemp_acc;
    std::vector<double> corr_avg = results.get_func("corr");
    corr_avg.shrink_to_fit();
    std::vector<double> vort_corr_avg = results.get_func("vort_corr");
    vort_corr_avg.shrink_to_fit();
    double N_sites = corr_avg.size();
    for (int i = 0; i < corr_avg.size(); ++i){
        corr_avg[i] /= (N_sites*N_sites);
        vort_corr_avg[i] /= (N_sites*N_sites);
    }
    loc_avg = mean(results.get_vals("loc"));
    loc2_avg = mean(results.get_vals("loc2"));
    loc4_avg = mean(results.get_vals("loc4"));
    loc_bind_avg = loc4_avg/loc2_avg/loc2_avg;
    loc_bind_err = bootstrap(results.get_vals("loc2"), 500, "binder");

    xmag_avg = mean(results.get_vals("xmag"));
    xmag2_avg = mean(results.get_vals("xmag2"));
    xmag4_avg = mean(results.get_vals("xmag4"));
    xmag_bind_avg = xmag4_avg/xmag2_avg/xmag2_avg;
    xmag_bind_err = bootstrap(results.get_vals("xmag2"), 500, "binder");

    mag_avg = mean(results.get_vals("mag"));
    mag2_avg = mean(results.get_vals("mag2"));
    mag4_avg = mean(results.get_vals("mag4"));
    mag_bind_avg = mag4_avg/mag2_avg/mag2_avg;
    mag_bind_err = bootstrap(results.get_vals("mag2"), 500, "binder");

    action_avg = mean(results.get_vals("action"));
    sx_avg = mean(results.get_vals("sx"));
    xkinks_avg = mean(results.get_vals("xkinks"));
    stagger_mag_avg = mean(results.get_vals("stagger_mag"));
    cluster_avg = mean(results.get_vals("cluster"));
    ptemp_acc = ptemp_moves/ptemp_total;

    std::cout << "Results for process #" << id << ":\n";
    timer.flag_end_time("full simulation");
    timer.print_timers();
    wolff.print_timers();
    std::cout << "total steps = " << num_step << ", acceptance ratio = " << ((double)num_accept)/((double)num_step)
            << ", number of measurements = " << num_meas
            << ", parallel tempering steps (or dump steps for no gpu) = " << ptemp_total
            << ", parallel tempering move probability: " << ptemp_acc
            << ", cluster size: " << cluster_avg/params.lengths[0]/params.lengths[1] << "\n\n\n";


//
//  Write dump files for the following items:
//      1. interactions matrix ("interactions#.csv")
//      2. timer info ("timers.csv")
//      3. list of measurements for all threads, including param info at the top ("dump#.csv")
//      4. final results table ("results.csv") - this can be done later after results are shared via MPI
//      5. correlation functions ("correlation#.csv")
//
	char dumpfile_name[100], corrfile_name[100], vortcorrfile_name[100], spacecorrfile_name[100], statefile_name[100];
	sprintf(dumpfile_name, "dump/dump%d.csv", id);
	sprintf(corrfile_name, "dump/correlation%d.csv", id);
    sprintf(vortcorrfile_name, "dump/vort_corr%d.csv", id);
    sprintf(statefile_name, "dump/state%d.pbm", id);
    sprintf(spacecorrfile_name, "dump/space_corr%d.csv", id);

	file.open(dumpfile_name);
    for (std::string name : results.names) {
        file << name << "," << vec2str(results.get_vals(name)) << "\n";
    }
    file.close();

    file.open(corrfile_name);
    file << lat.get_Lx() << "," << lat.get_Ly() << "\n";
    file << vec2str(corr_avg) << "\n";
    file.close();

    file.open(vortcorrfile_name);
    file << lat.get_Lx() << "," << lat.get_Ly() << "\n";
    file << vec2str(vort_corr_avg) << "\n";
    file.close();

    file.open(statefile_name);
    file << lat.get_Lx() << " " << lat.get_Ly() << "\n";
    file << lat.to_bitmap();
    file.close();

    file.open(spacecorrfile_name);
    results.write_results(&file, "space_corr");
    file.close();

//
//  Send results to master process for easy writing to a file
//
    if(id == 0){
        std::vector<double> ind_vars(p);//vector of independent variables (can change from alpha if this changes in earlier parts of the code)
        std::vector<double> locs(p), loc2s(p), loc4s(p), loc_binds(p), loc_bind_errs(p), 
                            mags(p), mag2s(p), mag4s(p), mag_binds(p), mag_bind_errs(p),
                            xmags(p), xmag2s(p), xmag4s(p), xmag_binds(p), xmag_bind_errs(p),
                            actions(p), sxs(p), xkinkss(p), stagger_mags(p), clusters(p), ptemp_ratios(p), traversals(p);
        ind_vars[0] = ind_var;

        locs[0] = loc_avg;
        loc2s[0] = loc2_avg;
        loc4s[0] = loc4_avg;
        loc_binds[0] = loc_bind_avg;
        loc_bind_errs[0] = loc_bind_err;

        mags[0] = mag_avg;
        mag2s[0] = mag2_avg;
        mag4s[0] = mag4_avg;
        mag_binds[0] = mag_bind_avg;
        mag_bind_errs[0] = mag_bind_err;

        xmags[0] = xmag_avg;
        xmag2s[0] = xmag2_avg;
        xmag4s[0] = mag4_avg;
        xmag_binds[0] = xmag_bind_avg;
        xmag_bind_errs[0] = xmag_bind_err;

        actions[0] = action_avg;
        sxs[0] = sx_avg;
        xkinkss[0] = xkinks_avg;
        stagger_mags[0] = stagger_mag_avg;
        clusters[0] = cluster_avg;
	ptemp_ratios[0] = ptemp_acc;
        traversals[0] = traversal;

        MPI_Status Stat;
        for (int i = 1; i < p; ++i){
            MPI_Recv(&(ind_vars[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(locs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc2s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc4s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc_binds[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(loc_bind_errs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(mags[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag2s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag4s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag_binds[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(mag_bind_errs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(xmags[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(xmag2s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(xmag4s[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(xmag_binds[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(xmag_bind_errs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

            MPI_Recv(&(actions[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(sxs[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(xkinkss[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(stagger_mags[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(clusters[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(ptemp_ratios[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);
            MPI_Recv(&(traversals[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat);

        }
        //cout << "binder cumulants: " << vec2str(loc_binds) << "\n";
        //write results
        FILE * file;
        file = fopen("results.dat", "w");
	    fprintf(file, "%-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s %-10.10s\n", 
                "alpha", "loc2", "xmag2", "mag2", "loc_bind", "xmag_bind", "mag_bind", "sx", "xkinks", "xstagger", "action", "cluster", "ptemp_acc", "traversal");
        for (int i = 0; i < p; ++i){
            fprintf(file, "%-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n", 
                    ind_vars[i], loc2s[i], xmag2s[i], mag2s[i], loc_binds[i], xmag_binds[i], mag_binds[i], sxs[i], xkinkss[i], stagger_mags[i], 
                    actions[i]/params.lengths[0]/params.lengths[1], clusters[i]/params.lengths[0]/params.lengths[1], ptemp_ratios[i], traversals[i]);
        }
        fclose(file);

        std::ofstream file2;

        file2.open("dump/results.csv");
        file2 << "ind variable," << vec2str(ind_vars) << "\n";
        file2 << "loc," << vec2str(locs) << "\n";
        file2 << "loc2," << vec2str(loc2s) << "\n";
        file2 << "loc4," << vec2str(loc4s) << "\n";
        file2 << "binder_y," << vec2str(loc_binds) << "\n";
        file2 << "binder_y_err," << vec2str(loc_bind_errs) << "\n";
        file2 << "xmag," << vec2str(xmags) << "\n";
        file2 << "xmag2," << vec2str(xmag2s) << "\n";
        file2 << "xmag4," << vec2str(xmag4s) << "\n";
        file2 << "binder_x," << vec2str(xmag_binds) << "\n";
        file2 << "binder_x_err," << vec2str(xmag_bind_errs) << "\n";
        file2 << "mag," << vec2str(mags) << "\n";
        file2 << "mag2," << vec2str(mag2s) << "\n";
        file2 << "mag4," << vec2str(mag4s) << "\n";
        file2 << "binder_tot," << vec2str(mag_binds) << "\n";
        file2 << "binder_tot_err," << vec2str(mag_bind_errs) << "\n";

        file2 << "sx," << vec2str(sxs) << "\n";
        file2 << "xkinks," << vec2str(xkinkss) << "\n";
        file2 << "stagger_mag," << vec2str(stagger_mags) << "\n";
        file2 << "action," << vec2str(actions) << "\n";
        file2 << "cluster," << vec2str(clusters) << "\n";
	file2 << "ptemp_acceptance," << vec2str(ptemp_ratios) << "\n";
        file2 << "traversal," << vec2str(traversals) << "\n";
        file2.close();
    }
    else{
        MPI_Send(&ind_var, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&loc_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc2_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc4_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc_bind_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&loc_bind_err, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&mag_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag2_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag4_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag_bind_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&mag_bind_err, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&xmag_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&xmag2_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&xmag4_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&xmag_bind_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&xmag_bind_err, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);

        MPI_Send(&action_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&sx_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&xkinks_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&stagger_mag_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&cluster_avg, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&ptemp_acc, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&traversal, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        
    }
	
//
//  End the MPI process
//
    if(mpi_use){
        MPI_Finalize();
    }

return 0;
}
