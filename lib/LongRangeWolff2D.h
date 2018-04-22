#pragma once
#include <thrust/host_vector.h>
#include <sstream>
#include "IsingLattice2D.h"
#include "class_mc_io.h"
#include "Matrix.h"
#include "MemTimeTester.h"
extern "C" {
#include "random.h"
}
#include <gsl/gsl_sf_zeta.h>
#include <cmath>
#include "obs_calc_fast.cuh"

double real_trigamma_cmplxarg(double x, double y, int m, int l);

class GeneralLRW{
	//This is a more lightweight version of the long range wolff class, designed to be
	//more general - hopefully to admit the use of user-defined "test_spins" functions
	//so that it can be used for many different MC simulations.  It does not hold the
	//lattice in order to keep the lattice easily available for GPU computations
	std::vector<std::vector<double>> interactions;//interaction matrix
	std::vector<std::vector<double>> interaction_sum;//values from interaction sum integral
	std::vector<double> site_sums;//absolute sum of interactions at each site
	double mag;//current magnetization
	double h; //applied field
	std::vector<spin> buffer;//holds the indices of spins to check
	std::vector<std::vector<int>> cluster;//0 at indices that are not in cluster, 1 at indices that are
	int cluster_size;//size of the cluster during a given MC step
	int cluster_mag;//total magnetization of the cluster
	void set_spin_boson_model_wc_exp(class_mc_params);//version for finite cutoff frequency with exponential tail
	void set_spin_boson_model_wc_hard(class_mc_params);//version for finite cutoff frequency with sharp cutoff
	void set_spin_boson_model(class_mc_params);//create interaction matrix for the spin boson model with infinite cutoff frequency
	void set_interaction_sum();
	void set_site_sums();
	MemTimeTester timer;

	inline double kernel_integrand(double alpha, int Nt, double rbv, double taub, double u){
		//evaluate the integrand of the kernel of the spin-boson effective interactions
		return 2.0 * alpha / Nt / Nt * u * cos(rbv*u)*(exp(u*(-taub)) + exp(u*(taub - 1.0)))/(1 - exp(-u));
	}

	double simpson_3_8(double alpha, int Nt, double bomc, double rbv, double taub, int N_points){
		//perform an integration of the spin-boson kernel with a hard, finite cutoff bomc.  Evaluated at 
		//site distance rbv = r/(beta*v), taub = tau/beta, with number of time slices Nt and number of 
		//integration slices N_points (this method will multiply N_points by 3 and use simpson's 3/8 rule)
		double du = (bomc)/(3.0*N_points);
		double result = 3.0*du*(4.0*alpha / Nt / Nt + kernel_integrand(alpha, Nt, rbv, taub, bomc))/8.0;
		for (int i = 1; i < 3*N_points; ++i){
			if(3*(i/3) == i){
				result += 0.75*du*kernel_integrand(alpha, Nt, rbv, taub, i*du);
			}
			else{
				result += 9.0*du*kernel_integrand(alpha, Nt, rbv, taub, i*du)/8.0;				
			}
		}
		return result;
	}

public:
	GeneralLRW(class_mc_params);

	void print_timers(){ timer.print_timers(); }

	void step(IsingLattice2D& lat);

    bool step_one_site(IsingLattice2D& lat, double *prev_action, thrust::host_vector<double>& corr_ref);

	bool step_one_site_prealloc(IsingLattice2D& lat, double *prev_action, thrust::host_vector<double>& corr_ref, 
			cufftHandle *full_forward_plan, cufftHandle *full_backward_plan, cufftHandle *onesite_forward_plan, cufftHandle *onesite_backward_plan,
			cudaError cuda_status, cufftDoubleReal *full_state_rs, cufftDoubleComplex *full_state_ft, cufftDoubleReal *onesite_state_rs, cufftDoubleComplex *onesite_state_ft);

	bool step_one_site_prealloc_lowmem(IsingLattice2D& lat, double *prev_action, thrust::host_vector<double>& corr_ref, 
			thrust::host_vector<double>& ,thrust::host_vector<double>& ,thrust::host_vector<double>& ,thrust::host_vector<double>& ,
			thrust::host_vector<double>& ,thrust::host_vector<double>& ,
			cufftHandle *full_forward_plan, cufftHandle *full_backward_plan, cufftHandle *onesite_forward_plan, cufftHandle *onesite_backward_plan,
			cudaError cuda_status, cufftDoubleReal *full_state_rs, cufftDoubleComplex *full_state_ft, cufftDoubleReal *onesite_state_rs, cufftDoubleComplex *onesite_state_ft);

	bool metropolis_step(IsingLattice2D& lat, double *action);

	bool fast_metropolis_step(IsingLattice2D& lat, double *action);

	void test_spins(spin, IsingLattice2D& lat);

    void test_spins_one_site(spin, IsingLattice2D& lat);

	void set_mag(IsingLattice2D& lat){
		mag = 0;
		for (int i = 0; i < lat.get_Lx(); ++i){
			for (int j = 0; j < lat.get_Ly(); ++j){
				mag += lat.get_spin(i, j);
			}
		}
		mag = ((double) mag) / lat.get_N();
	}

	double get_mag(){
		return mag;
	}

	double get_cluster_size() { return (double)cluster_size; }

    void get_thrust_interactions(thrust::host_vector<double>& result){
		if(result.size() != interactions.size()*interactions[0].size()){
			result.resize(interactions.size()*interactions[0].size());
		}
		for (int i = 0; i < interactions.size(); ++i){
			for (int j = 0; j < interactions[i].size(); ++j){
				result[i*interactions[i].size() + j] = interactions[i][j];
			}
		}
    }

	std::vector<double> get_interaction_vector(){
		std::vector<double> result(interactions.size()*interactions[0].size());
		int Ly = interactions[0].size();
		for (int x = 0; x < interactions.size(); ++x){
			for (int y = 0; y < Ly; ++y){
				result[x*Ly + y] = interactions[x][y];
			}
		}
		return result;
	}

	std::vector<std::vector<double>> get_interactions(){
		return interactions;
	}
    
	double calc_sx(IsingLattice2D&);

	double calc_space_kinks(IsingLattice2D& lat);

	double calc_sz_stagger(IsingLattice2D&);

	double calc_xmag(IsingLattice2D&);

    double calc_mag(IsingLattice2D&);

	double calc_loc(IsingLattice2D&);

	double calc_loc2(IsingLattice2D&);

	double calc_loc4(IsingLattice2D&);

	double calc_s1s2(IsingLattice2D&);

	double calc_action_slow(IsingLattice2D&);

	double calc_point_action_slow(IsingLattice2D&,int, int);

	double calc_point_action_fast(IsingLattice2D&, int, int);

	void print_site_sums(){
		std::cout << "Site sums:\n" << vec2str(site_sums) << "\n";
	}

	void print_interactions() {
		std::stringstream outstring;
		outstring << "Interactions:\n";
		for (int i = 0; i < interactions.size(); ++i){
			outstring << vec2str(interactions[i]) << "\n";
		}
		std::cout << outstring.str() << "\n";
	}


	void print_interaction_sum() {
		std::stringstream outstring;
		outstring << "Interaction sum:\n";
		for (int i = 0; i < interaction_sum.size(); ++i){
			outstring << vec2str(interaction_sum[i]) << "\n";
		}
		std::cout << outstring.str() << "\n";
	}

	std::string get_int_string() {
		std::stringstream ss;
		ss << "Interactions:\n" ;
		for (int i = 0; i < interactions.size(); ++i){
			ss << vec2str(interactions[i]) << "\n";
		}
		return ss.str();
	}

};
