#include "MemTimeTester.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "obs_calc_fast.cuh"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/count.h>
#include <thrust/inner_product.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <string>
#include <sstream>


#include <stdio.h>


int main(int argc, char* argv[])
{
	//first input must be the number of states in the dump directory
	std::stringstream input(argv[1]);
	int num_states;
	input >> num_states;

	MemTimeTester timer;
	std::ifstream file;
	std::string line;
	int Lx, Ly;


	//check average energy, correlation function
	std::cout << "Calculating average of " << num_states << " runs energy and correlation\n\n";
	char filename[300];

	thrust::host_vector<int> state;
	thrust::host_vector<int>& state_ref = state;
	thrust::host_vector<double> corr;
	thrust::host_vector<double>& corr_ref = corr;
	thrust::device_vector<double> corr_temp = corr;
	thrust::device_vector<double> corr_total;
	for (int run = 0; run < num_states; ++run) {
		sprintf(filename, "./dump/state%d.csv", run);
		file.open(filename);
		if (file.is_open()) {
			file >> Lx >> Ly;
			if(run == 0){
				state.resize(Lx*Ly);
			}
			for (int i = 0; i < Lx; ++i) {
				for (int j = 0; j < Ly; ++j) {
					file >> state[i*Ly + j];
				}
			}
		}
		else {
			std::cout << "Error: input file not opened\n";
		}
		file.close();
		if (run == 0){
			corr.resize(Lx*Ly);
			corr_temp.resize(Lx*Ly);
			corr_total.resize(Lx*Ly);
			thrust::fill(corr_total.begin(), corr_total.end(), 0.0);
		}
		timer.flag_start_time("double site correlation measurement");
		if(Lx == 1){
			calc_corr_fast_1site(corr_ref, state_ref, Ly);
		}
		else {
			calc_corr_fast_2site(corr_ref, state_ref, Ly);
		}
		//corr_temp = corr_ref;
		thrust::copy(corr_ref.begin(), corr_ref.end(), corr_temp.begin());
		thrust::transform(corr_temp.begin(), corr_temp.end(), corr_total.begin(), corr_total.begin(), thrust::plus<double>());
		timer.flag_end_time("double site correlation measurement");
	}
	thrust::constant_iterator<double> factor(1.0 / num_states);
	thrust::transform(corr_total.begin(), corr_total.end(), factor, corr_temp.begin(), thrust::multiplies<double>());
	//corr = corr_temp;
	thrust::copy(corr_temp.begin(), corr_temp.end(), corr.begin());
/*
	std::cout << "Double Site Correlation function:\n";
	for (int i = 0; i < corr.size(); ++i) {
		std::cout << corr[i] << ",";
	}
	std::cout << "\n";
*/
	std::ofstream outfile;
	outfile.open("corr.csv");
	for (int i = 0; i < Lx; ++i) {
		for (int j = 0 ; j < Ly - 1; ++j){
               		outfile  << corr[i*Ly + j] << ",";
		}
		outfile << corr[i * Ly + Ly - 1] << "\n";
        }
	outfile.close();

	timer.print_timers();

    return 0;
}
