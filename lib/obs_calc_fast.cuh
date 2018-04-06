#pragma once
#include <fstream>
#include <vector>
#include <thrust/host_vector.h>
#include "cufft.h"
#include "MemTimeTester.h"
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/copy.h>
#include <thrust/reduce.h>
#include <thrust/gather.h>
#include <thrust/scan.h>
#include <thrust/fill.h>

typedef float2 Complex;

void show_memory();

template <typename InputIterator1,
	typename InputIterator2,
	typename OutputIterator>
	OutputIterator expand(InputIterator1 first1,
		InputIterator1 last1,
		InputIterator2 first2,
		OutputIterator output);

void calc_corr_fast_1site(thrust::host_vector<double>&, thrust::host_vector<int>&, int);

void calc_corr_fast_2site(thrust::host_vector<double>&, thrust::host_vector<int>&, int);

double calc_action_fast(thrust::host_vector<double>& corr, thrust::host_vector<double>& interactions);

__global__ void calc_shift_state(double* state, double* shift_state, int i, int j);

double thrust_calc_action_general(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, int Lx, int Ly, dim3 threads);

__global__ void elementwise_product_cmplx(cufftDoubleComplex *source);

double cufft_calc_action(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, thrust::host_vector<double>& corr, int Lx, int Ly);

double cufft_calc_action_timer(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, thrust::host_vector<double>& corr, int Lx, int Ly, MemTimeTester * timer);

double cufft_calc_action_timer_prealloc(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions,
	thrust::host_vector<double>& corr, int Lx, int Ly, MemTimeTester * timer,
	cufftHandle forward_plan, cufftHandle backward_plan, cudaError cuda_status, cufftDoubleReal *state_rs, cufftDoubleComplex *state_ft);

double cufft_calc_action_prealloc(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions,
	thrust::host_vector<double>& corr, int Lx, int Ly,	cufftHandle forward_plan, cufftHandle backward_plan, cudaError cuda_status, cufftDoubleReal *state_rs, cufftDoubleComplex *state_ft);

__global__ void point_action_shift(double *state, double *shift_state, int x, int y, int Lx, int Ly);

double calc_point_action(thrust::host_vector<double>& state, thrust::host_vector<double>& interactions, int x, int y, int Lx, int Ly);



