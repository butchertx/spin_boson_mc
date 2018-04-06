#include <iostream>
#include <fstream>
#include "MemTimeTester.h"
#include "IsingLattice2D.h"
#include "LongRangeWolff2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
extern "C" {
#include "random.h"
}

int main(int argc, char* argv[]) {

	//Get input parameters
	class_mc_params params;
	std::ifstream infile;
	std::cout << "Program name: " << argv[0] << "\n";
	std::cout << "Input file: " << argv[1] << "\n\n";
	infile.open(argv[1]);
	read_input_ising(&infile, &params);
	read_input_spin_boson(&infile, &(params.sbparams));
	apply_spin_boson_params(&params);
	std::cout << "Parameters: \n" << params.to_string();
	write_params(params);


	//Set up lattice, wolff, and measurement objects
	IsingLattice2D test_lat = IsingLattice2D(params.lengths[0], params.lengths[1], params.rand_seed);
	LongRangeWolff2D wolff = LongRangeWolff2D(&test_lat, &params);
	class_mc_measurements results;
	MemTimeTester timer;
	results.names = { "mag", "mag2", "mag4", "cluster size", "a_mag", "a_mag2", "sx", "loc"};
	results.values = { {}, {}, {} , {}, {}, {}, {}, {}};
	
	//Set flags within wolff algorithm
	if (params.alg.compare("long_range_cluster") == 0) { wolff.set_alg_long_range_cluster(); }
	else if (params.alg.compare("nearest_neighbor_cluster") == 0) { wolff.set_alg_short_range_cluster(); }
	wolff.set_alg_long_range_cluster();
	if(params.sbparams.A0 == 0){
		wolff.set_alg_short_range_cluster();
	}
	wolff.write_interactions();
	
	//Set up parallel tempering lattices
	int num_ptemp = 1;
	std::cout << "Number of parallel tempering lattices: " << num_ptemp << "\n";
	std::vector<LongRangeWolff2D> wolffs = { wolff };


	//Run simulation
	timer.flag_start_time("simulation");
		//equilibrate
	for (int i = 0; i < params.eq_time; ++i) {
		for (int ptemp = 0; ptemp < num_ptemp; ++ptemp) {
			wolffs[ptemp].step();
		}
	}
		//measure
	double mag, a_mag;
	IsingLattice2D lat_copy = wolffs[0].get_lat();
	IsingLattice2D& lat_ref = lat_copy;
	for (int dump = 0; dump < params.max_dumps; ++dump) {
		for (int measure = 0; measure < params.measures_per_dump; ++measure) {
			for (int step = 0; step < params.steps_per_measure; ++step) {
				for (int ptemp = 0; ptemp < num_ptemp; ++ptemp) {
					wolffs[ptemp].step();
				}
			}
			//make measurements
			mag = wolffs[0].get_mag();
			a_mag = wolffs[0].calc_sz_stagger();
/*
			if (mag > 1.0){
				wolffs[0].output_state();
				std::cout << "recorded mag: " << mag << "\n";
			}
*/
			results.record("mag", mag);
			results.record("mag2", mag * mag);
			results.record("mag4", mag*mag*mag*mag);
			results.record("cluster size", wolffs[0].get_cluster_size());
			results.record("a_mag", a_mag);
			results.record("a_mag2", a_mag*a_mag);
			results.record("sx", wolffs[0].calc_sx());
			results.record("loc", wolffs[0].calc_localization());
/*
			timer.flag_start_time("state_write");
			lat_copy = wolffs[0].get_lat();
			write_state(dump*params.measures_per_dump + measure, lat_ref);
			timer.flag_end_time("state_write");
*/
		}
		write_outputs_var(dump, results);
		std::cout << "Dump " << dump + 1 << " out of " << params.max_dumps << "\n";
	}
	write_final_outputs(results, params.max_dumps);
	timer.flag_end_time("simulation");

	timer.print_timers();
}
