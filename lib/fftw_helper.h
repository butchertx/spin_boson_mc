#pragma once
#include <fftw3.h>
#include <vector>

void dagger_norm(fftw_complex* s_k, fftw_complex* c_k, int N);

std::vector<double> calc_correlation_no_wisdom(std::vector<double> state, int Lx, int Ly);

std::vector<double> calc_skk_no_wisdom(std::vector<int>& state);

double calc_action_fft(std::vector<double>& corr, std::vector<double> interactions, std::vector<double> state, int Lx, int Ly);

double calc_action_slow(std::vector<double>& interactions, std::vector<double> state);

