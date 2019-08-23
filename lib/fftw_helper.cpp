#include "fftw_helper.h"
#include <vector>
#include <iostream>

#ifdef TEST
    int main(){
        std::cout << "Testing Helper functions\n";
        std::vector<double> state1 = {1,1,1,1,1,1,1,1,1,1};
        std::vector<double> state2 = {1,-1,1,-1,1,-1,1,-1,1,-1};
        std::vector<double> interactions = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        std::vector<double> corr(10);
        double slow_action1 = calc_action_slow(interactions, state1);
        double slow_action2 = calc_action_slow(interactions, state2);
        double fast_action1 = calc_action_fft(corr, interactions, state1, 2, 5);
        double fast_action2 = calc_action_fft(corr, interactions, state2, 2, 5);

        std::cout << "Results:\n"
            << "Slow Action Ferro = " << slow_action1 << "\n"
            << "Fast Action Ferro = " << fast_action1 << "\n"
            << "Slow Action Antiferro = " << slow_action2 << "\n"
            << "Fast Action Antiferro = " << fast_action2 << "\n";
        return 0;
    }
#endif

void dagger_norm(fftw_complex* s_k, fftw_complex* c_k, int N){
    //c_k = s(k) dot s(-k) = s(k) dot s*(k)
    for(int i = 0; i < N; ++i){
        c_k[i][0] = s_k[i][0]*s_k[i][0] + s_k[i][1]*s_k[i][1];
        c_k[i][1] = 0;
    }
}

std::vector<double> calc_correlation_no_wisdom(std::vector<double> state, int Lx, int Ly){
    int L = state.size();
    std::vector<double> result(L);
    fftw_complex *in, *out;
    fftw_plan p_forward, p_backward;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    for(int i = 0; i < L; ++i){
        in[i][0] = state[i];
        in[i][1] = 0.0;
    }
    p_forward = fftw_plan_dft_2d(Lx, Ly, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p_backward = fftw_plan_dft_2d(Lx, Ly, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);
    dagger_norm(out, in, L);
    fftw_execute(p_backward);
    for(int i = 0; i < L; ++i){
        result[i] = out[i][0]/L;
    }
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(in); fftw_free(out);

    return result;
}

std::vector<double> calc_skk_no_wisdom(std::vector<int>& state){
    int L = state.size();
    std::vector<double> result(L);
    fftw_complex *in, *out;
    fftw_plan p_forward;//, p_backward;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);
    for(int i = 0; i < L; ++i){
        in[i][0] = state[i];
    }
    p_forward = fftw_plan_dft_1d(L, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    //p_backward = fftw_plan_dft_1d(L, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);
    dagger_norm(out, in, L);
    //fftw_execute(p_backward);
    for(int i = 0; i < L; ++i){
        result[i] = out[i][0]/L;
    }
    fftw_destroy_plan(p_forward);
    //fftw_destroy_plan(p_backward);
    fftw_free(in); fftw_free(out);

    return result;
}

double calc_action_fft(std::vector<double>& corr, std::vector<double> interactions, std::vector<double> state, int Lx, int Ly){
    corr = calc_correlation_no_wisdom(state, Lx, Ly);
    double result = 0.0;
    for(int i = 0; i < state.size(); ++i){
        result += 0.5*interactions[i]*corr[i];
    }
    return result;
}

double calc_action_slow(std::vector<double>& interactions, std::vector<double> state){
    double result = 0;
    int L = state.size();
    for(int i = 0; i < L; ++i){
        for(int x = 1; x < L; ++x){
            result += 0.5*state[i]*state[(i+x)%L]*interactions[x];
        }
    }
    return result;
}
