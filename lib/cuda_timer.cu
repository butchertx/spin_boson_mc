
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <thrust/host_vector.h>
#include "MemTimeTester.h"
#include "obs_calc_fast.cuh"
extern "C" {
#include "random.h"
}
/**
For my GPU, the maximum threads per multiprocessor is 2048, max per block is 1024, and the max dimensions are (1024, 1024, 64)
Maximum shared memory per block is 49152 bytes, L2 Cache is 2097152 bytes, 15 multiprocessors, 128 CUDA cores/MP
Cuda capability: 6.1, CUDA driver version 8.0
**/

double calc_action_slow(thrust::host_vector<double>& lat, thrust::host_vector<double>& interactions, int Lx, int Ly) {
	//essentially, the time-averaged energy of a given set of quantum fluctuations.  Calculate the action then divide by beta
	double S = 0;
	double s1;
	for (int i = 0; i < Lx; ++i) {
		for (int j = 0; j < Ly; ++j) {
			s1 = lat[i*Ly + j];
			for (int m = 0; m < Lx; ++m) {
				for (int n = 0; n < Ly; ++n) {
					S += s1 * lat[m*Ly + n] * interactions[((m - i + Lx) % Lx)*Ly + ((n - j + Ly) % Ly)];
				}
			}
		}
	}

	return 0.5*S;
}

thrust::host_vector<double> rand_vector(int length) {
	thrust::host_vector<double> result(length);
	for (int i = 0; i < length; ++i) {
		result[i] = drand1_();
	}
	return result;
}

thrust::host_vector<double> transpose(thrust::host_vector<double> s, int Lx, int Ly) {
	//take a vector indexed like s(x, y) = s[x*Ly + y] and make it s(x,y) = s[y*Lx + x]
	thrust::host_vector<double> new_s(Lx*Ly);
	for (int x = 0; x < Lx; ++x) {
		for (int y = 0; y < Ly; ++y) {
			new_s[y*Lx + x] = s[x*Ly + y];
		}
	}
	return new_s;
}

int main() {
	std::cout << "Timing different versions of fast calc and checking for accuracy of results\n";
	MemTimeTester timer;
	int seed = 1892347;
	rand_init_(&seed);
	double fast_action;
	double slow_action;
	
	/*Test fast action calculations
	*	1. Test accuracy for small calculations
	*	2. Test speed for some different versions
	*/

	//1. Test accuracy
	thrust::host_vector<double> acc_state = rand_vector(128);//just do something simple like lx = 8, ly = 16
	thrust::host_vector<double>& acc_state_ref = acc_state;
	thrust::host_vector<double> acc_int = rand_vector(128);//same dimensions
	thrust::host_vector<double>& acc_int_ref = acc_int;
	thrust::host_vector<double> acc_corr(128, 0.0);
	thrust::host_vector<double>& acc_corr_ref = acc_corr;
	slow_action = calc_action_slow(acc_state_ref, acc_int_ref, 8, 16);
	fast_action = cufft_calc_action(acc_state_ref, acc_int_ref, acc_corr_ref, 8, 16);
	if (slow_action == fast_action) {
		std::cout << "Action accuracy test passed!\n\n";
	}
	else {
		std::cout << "Action accuracy not quite right! Slow action: " 
			<< slow_action << ", fast action: " << fast_action << ", abs difference: " << abs(fast_action - slow_action) << "\n\n";
	}

	//2. Test timing
	thrust::host_vector<double> one_state(32 * 1024, 1.0);
	thrust::host_vector<double> one_state_2(32 * 1024, 1.0);
	thrust::host_vector<double> corr(32 * 1024, 1.0);
	thrust::host_vector<double>& corr_ref = corr;
	thrust::host_vector<double>& one_int_ref = one_state_2;
	thrust::host_vector<double>& one_state_ref = one_state;
	int Lx = 32, Ly = 1024;
	//no preallocation
	timer.flag_start_time("cufft total no prealloc");
	for (int i = 0; i < 10000; ++i) {
		fast_action = cufft_calc_action_timer(one_state_ref, one_int_ref, corr_ref, Lx, Ly, &timer);
		cudaThreadSynchronize();
		timer.flag_end_time("thrust inner product");
	}
	timer.flag_end_time("cufft total no prealloc");
	if (fast_action == 0.5 * 32 * 32 * 1024 * 1024) {
		std::cout << "Fast no prealloc action correct!\n";
	}
	else {
		std::cout << "Fast no prealloc action: " << fast_action << ", should be: " << 0.5 * 32 * 32 * 1024 * 1024 << "\n";
	}



	//with preallocation
	timer.flag_start_time("cufft total with prealloc");
	cufftHandle forward_plan, backward_plan;
	cudaError cuda_status;
	cufftDoubleReal *state_rs;//real space state
	cufftDoubleComplex *state_ft;//fourier space state
	int n[2] = { Lx, Ly };
	//allocate real state
	cuda_status = cudaMalloc((void**)&state_rs, sizeof(cufftDoubleReal)*Lx*Ly);
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to allocate. Code: %s\n", cudaGetErrorName(cuda_status));
		show_memory();
		return 0;
	}

	//allocate transform state
	cuda_status = cudaMalloc((void**)&state_ft, sizeof(cufftDoubleComplex)*Lx*(Ly / 2 + 1));
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Cuda error: Failed to allocate. Code: %s\n", cudaGetErrorName(cuda_status));
		show_memory();
		return 0;
	}
	cudaThreadSynchronize();
	timer.flag_end_time("cufft allocation");

	timer.flag_start_time("cufft prealloc plan");
	//Create 2D R2C FFT plan
	if (cufftPlanMany(&forward_plan, 2, n, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, 1) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT Error: Unable to create plan\n");
		return 0;
	}
	cudaThreadSynchronize();
	timer.flag_end_time("cufft prealloc plan");


	timer.flag_start_time("cufft prealloc plan");
	//Inverse fourier transform plan
	if (cufftPlanMany(&backward_plan, 2, n, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, 1) != CUFFT_SUCCESS) {
		fprintf(stderr, "CUFFT Error: Unable to create plan\n");
		return 0;
	}
	cudaThreadSynchronize();
	timer.flag_end_time("cufft prealloc plan");

	for (int i = 0; i < 10000; ++i) {
		fast_action = cufft_calc_action_timer_prealloc(one_state_ref, one_int_ref, corr_ref, Lx, Ly, &timer,
			forward_plan,  backward_plan,  cuda_status,  state_rs,  state_ft);
	}

	timer.flag_start_time("cufft prealloc deallocate");
	//Deallocate memory
	cufftDestroy(forward_plan);
	cufftDestroy(backward_plan);
	cudaFree(state_rs);
	cudaFree(state_ft);
	cudaThreadSynchronize();
	timer.flag_end_time("cufft prealloc deallocate");
	timer.flag_end_time("cufft total with prealloc");

	if (fast_action == 0.5 * 32 * 32 * 1024 * 1024) {
		std::cout << "Fast prealloc action correct!\n";
	}
	else {
		std::cout << "Fast prealloc action: " << fast_action << ", should be: " << 0.5 * 32 * 32 * 1024 * 1024 << "\n";
	}

	timer.print_timers();

	return 0;
}