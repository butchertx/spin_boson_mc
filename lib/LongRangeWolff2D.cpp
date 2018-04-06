#define _USE_MATH_DEFINES
#include <cmath>
#include <omp.h>
#include <algorithm>
#include "LongRangeWolff2D.h"

double real_trigamma_cmplxarg(double x, double y, int m, int l){
    //compute the real part of the trigamma function at z = x + iy using
    //my accelerated series formula
    //argument should have positive real part (x > 0)
    //m gives the number of zeta function terms to use
    //l gives the number of terms to use in the residual series
    double result = 0;

    //if x is small, the result will be large and thus, inaccurate.  Use the
    //polygamma shift formula to give a larger value of x for the computation
    if(x < 1000){
        return ((x*x - y*y)/(x*x + y*y)/(x*x + y*y)) + real_trigamma_cmplxarg(x + 1, y, m, l);
    }
    else{
        double phase, ypow = 1.0, xpow, denom;
        //compute finite sum of Hurwitz zeta functions
        for (int i = 1; i < m; ++i){
            phase = 2*(i/2) == i ? -1.0 : 1.0;
            result += phase*(2*i - 1)*ypow*gsl_sf_hzeta(2*i, x);
            ypow = ypow*y*y;
        }
        //compute the infinite sum of residuals
        phase = 2*(m/2) == m ? 1.0 : -1.0;
        for (int n = 0; n < l; ++n){
            xpow = pow(x + n, 2*m);
            denom = ((x + n)*(x + n) + y*y)*((x + n)*(x + n) + y*y);
            result += phase*ypow*((2*m - 1)*y*y + (2*m + 1)*(x + n)*(x + n)) / xpow / denom;
        }
        return result;
    }
}

/**
----------------------------------------------------------------------------------------------------
General LRW
----------------------------------------------------------------------------------------------------
**/

GeneralLRW::GeneralLRW(class_mc_params params_in) {
	//generate member objects
	timer.flag_start_time("Setup");
	for(int i = 0; i < params_in.lengths[0]; ++i){
		interactions.push_back({});
		interaction_sum.push_back({});
		cluster.push_back({});
		for (int j = 0; j < params_in.lengths[1]; ++j){
			interactions[i].push_back(0);
			interaction_sum[i].push_back(0);
			cluster[i].push_back(0);
		}
	}
	cluster_size = 0;
	cluster_mag = 0;

	if (params_in.cutoff_type.compare("exp") == 0){
		//std::cout << "Using Finite Exponential Cutoff\n";
		set_spin_boson_model_wc_exp(params_in);
	}
	else if (params_in.cutoff_type.compare("inf") == 0){
		//std::cout << "Using Infinite Cutoff\n";
		set_spin_boson_model(params_in);
	}
	else{
		//std::cout << "Using Finite Hard Cutoff\n";
		set_spin_boson_model_wc_hard(params_in);
	}
	set_interaction_sum();
	set_site_sums();
	mag = params_in.lengths[0]*params_in.lengths[1] + 1;
	h = params_in.h;
	timer.flag_end_time("Setup");
}

void GeneralLRW::set_spin_boson_model(class_mc_params params){
	//Infinite cutoff
	//model is translationally invariant, so (i,j) will represent a distance in the (spacial, imaginary time) direction
	double g = params.sbparams.g, A = params.sbparams.A0, gamma = params.Js[1], J = params.Js[0], v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], x, y;
	int Nt = params.lengths[1], Lx = params.lengths[0];
	for (int i = 0; i < interactions.size(); ++i) {
		for (int j = 0; j < interactions[i].size(); ++j) {
			if (i == 0) {
				//temporal nearest neighbor interaction
				if (j == 1 || j == Nt - 1) {
					interactions[i][j] = -gamma - 0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI / Nt) / sin(M_PI / Nt);
				}
				//temporal self-interaction
				if (j > 1 && j < Nt - 1) {
					interactions[i][j] = -0.5*A*M_PI*M_PI / Nt / Nt / sin(M_PI * j / Nt) / sin(M_PI * j / Nt);
				}
			}
			else if (j == 0) {
				//spacial nearest neighbor interactions
				if (i == 1 || i == Lx - 1) {
					interactions[i][j] = -J + A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a / Nt / v / tc) / sinh(M_PI * a / Nt / v / tc);
				}
				//spacial same-time long range interactions
				if (i > 1 && i < Lx - 1) {
					if (i > Lx / 2){
						interactions[i][j] = A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a * (Lx - i) / Nt / v / tc) / sinh(M_PI * a * (Lx - i) / Nt / v / tc);
					}
					else{
						interactions[i][j] = A*M_PI*M_PI / Nt / Nt / sinh(M_PI * a * i / Nt / v / tc) / sinh(M_PI * a * i / Nt / v / tc);
					}
				}
			}
			else {
				//different site, different time interactions
				x = M_PI*j / Nt;
				if (i > Lx / 2){
					y = M_PI*(Lx - i)*a / v / Nt / tc;;
				}
				else{
					y = M_PI*i*a / v / Nt / tc;
				}
				interactions[i][j] = -M_PI*M_PI*A / Nt / Nt*(sin(x)*sin(x)*cosh(y)*cosh(y) - cos(x)*cos(x)*sinh(y)*sinh(y)) /
					((sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y))*(sin(x)*sin(x)*cosh(y)*cosh(y) + cos(x)*cos(x)*sinh(y)*sinh(y)));
			}
		}
	}
}

void GeneralLRW::set_spin_boson_model_wc_exp(class_mc_params params){
	//model is translationally invariant, so (i,j) will represent a distance in the (spacial, imaginary time) direction
	double g = params.sbparams.g, A = params.sbparams.A0, gamma = params.Js[1], J = params.Js[0], v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], wc = params.sbparams.omega_c, x, y;
	int Nt = params.lengths[1], Lx = params.lengths[0];
	double beta = tc*Nt;
	double omcb = 1.0 / wc / beta;//1/beta/omega_c
	double prefactor = A / Nt / Nt;
	for (int i = 0; i < interactions.size(); ++i) {
		y = (i > Lx / 2) ? a * (Lx - i) / v / beta : y = a * i / v / beta;//periodic boundary condition in spatial direction
		for (int j = 0; j < interactions[i].size(); ++j) {
			x = ((double)j) / Nt;//scaled time coordinate
			if (i == 0) {
				//temporal nearest neighbor interaction
				if (j == 1 || j == Nt - 1) {
					interactions[i][j] = -gamma - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10));
				}
				//temporal self-interaction
				if (j > 1 && j < Nt - 1) {
					interactions[i][j] = - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10));
				}
			}
			else if (j == 0) {
				//spacial nearest neighbor interactions
				if (i == 1 || i == Lx - 1) {
					interactions[i][j] = - J - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10)
																+ real_trigamma_cmplxarg(x + omcb, -y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, -y, 10, 10));
				}
				//spacial same-time long range interactions
				if (i > 1 && i < Lx - 1) {
					interactions[i][j] = - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10)
																+ real_trigamma_cmplxarg(x + omcb, -y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, -y, 10, 10));
				}
			}
			else {
				//different site, different time interactions
				interactions[i][j] = - prefactor * (real_trigamma_cmplxarg(x + omcb, y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, y, 10, 10)
																+ real_trigamma_cmplxarg(x + omcb, -y, 10, 10) + real_trigamma_cmplxarg(1 - x + omcb, -y, 10, 10));
			}
		}
	}
}

void GeneralLRW::set_spin_boson_model_wc_hard(class_mc_params params){
	//model is translationally invariant, so (i,j) will represent a distance in the (spacial, imaginary time) direction
	double g = params.sbparams.g, A = params.sbparams.A0, gamma = params.Js[1], J = params.Js[0], v = params.sbparams.v, a = params.spacings[0], tc = params.spacings[1], wc = params.sbparams.omega_c, x, y;
	int Nt = params.lengths[1], Lx = params.lengths[0];
	double beta = tc*Nt;
	double omcb = 1.0 / wc / beta;//1/beta/omega_c
	double prefactor = A / Nt / Nt;
	for (int i = 0; i < interactions.size(); ++i) {
		y = (i > Lx / 2) ? a * (Lx - i) / v / beta : y = a * i / v / beta;//periodic boundary condition in spatial direction
		for (int j = 0; j < interactions[i].size(); ++j) {
			x = ((double)j) / Nt;//scaled time coordinate
			if (i == 0) {
				//temporal nearest neighbor interaction
				if (j == 1 || j == Nt - 1) {
					interactions[i][j] = -gamma - 0.5*simpson_3_8(A, Nt, 1.0/omcb, y, x, 10000);
				}
				//temporal self-interaction
				if (j > 1 && j < Nt - 1) {
					interactions[i][j] = - 0.5*simpson_3_8(A, Nt, 1.0/omcb, y, x, 10000);
				}
			}
			else if (j == 0) {
				//spacial nearest neighbor interactions
				if (i == 1 || i == Lx - 1) {
					interactions[i][j] = - J - simpson_3_8(A, Nt, 1.0/omcb, y, x, 10000);
				}
				//spacial same-time long range interactions
				if (i > 1 && i < Lx - 1) {
					interactions[i][j] = - simpson_3_8(A, Nt, 1.0/omcb, y, x, 10000);
				}
			}
			else {
				//different site, different time interactions
				interactions[i][j] = - simpson_3_8(A, Nt, 1.0/omcb, y, x, 10000);
			}
		}
	}
}

void GeneralLRW::set_interaction_sum(){
	double total = 0;
	for (int i = 0; i < interaction_sum.size(); ++i){
		for (int j = 0; j < interaction_sum[i].size(); ++j){
			total += interactions[i][j] > 0 ? interactions[i][j] : -interactions[i][j];
			interaction_sum[i][j] = total;
		}
	}
}

void GeneralLRW::set_site_sums(){
	//THIS ASSUMES INTERACTION SUM HAS ALREADY BEEN SET
	//CHECK THIS
	for (int i = 0; i < interaction_sum.size(); ++i){
		site_sums.push_back(interaction_sum[i][interaction_sum[i].size() - 1]);
	}
}

double GeneralLRW::calc_sx(IsingLattice2D& lat) {
	//assume a quantum-classical mapping.  Then, <flux>_i = 1/beta * sum_t{1/2(1 - <s(t + dt)s(t)>)}
	double final_result = 0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		double site_mag = 0;
		for (int j = 0; j < lat.get_Ly() - 1; ++j) {
			site_mag += 0.5*(1 - lat.get_spin(i, j)*lat.get_spin(i, j + 1));
		}
		site_mag += 0.5*(1 - lat.get_spin(i, 0)*lat.get_spin(i, lat.get_Ly() - 1));
		final_result += site_mag / lat.get_Ly();
	}
	return final_result / lat.get_Lx();
}

double GeneralLRW::calc_space_kinks(IsingLattice2D& lat) {
	double final_result = 0;
	for (int i = 0; i < lat.get_Ly(); ++i) {
		double time_mag = 0;
		for (int j = 0; j < lat.get_Lx() - 1; ++j) {
			time_mag += 0.5*(1 - lat.get_spin(j, i)*lat.get_spin(j + 1,i));
		}
		time_mag += 0.5*(1 - lat.get_spin(0, i)*lat.get_spin(lat.get_Lx() - 1, i));
		final_result += time_mag / lat.get_Lx();
	}
	return final_result / lat.get_Ly();
}

double GeneralLRW::calc_xmag(IsingLattice2D& lat) {
	//spatial analogue of localization
	double final_result = 0.0;
	for (int j = 0; j < lat.get_Ly(); ++j) {
		double time_mag = 0.0;
		for (int i = 0; i < lat.get_Lx(); ++i) {
			time_mag += lat.get_spin(i, j);
		}
		final_result += (time_mag > 0) ? time_mag / lat.get_Lx() : -time_mag / lat.get_Lx();
	}
	return final_result / lat.get_Ly();
}

double GeneralLRW::calc_loc(IsingLattice2D& lat) {
	//calculate the localization - absolute value of sz at each site, averaged
	double final_result = 0.0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		double site_mag = 0.0;
		for (int j = 0; j < lat.get_Ly(); ++j) {
			site_mag += lat.get_spin(i, j);
		}
		site_mag = site_mag > 0 ? site_mag : -site_mag;
		final_result += site_mag / lat.get_Ly();
	}
	return final_result / lat.get_Lx();
}

double GeneralLRW::calc_mag(IsingLattice2D& lat){
    double tempmag = 0;
    for (int i = 0; i < lat.get_Lx(); ++i){
        for (int j = 0; j < lat.get_Ly(); ++j){
            tempmag += lat.get_spin(i, j);
        }
    }
    return ((double) tempmag) / lat.get_N();
}

double GeneralLRW::calc_loc2(IsingLattice2D& lat) {
	//calculate the localization - absolute value of sz at each site, squared and averaged
	double final_result = 0.0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		double site_mag = 0.0;
		for (int j = 0; j < lat.get_Ly(); ++j) {
			site_mag += lat.get_spin(i, j);
		}
		site_mag *= site_mag / lat.get_Ly() / lat.get_Ly();
		final_result += site_mag;
	}
	return final_result / lat.get_Lx();
}

double GeneralLRW::calc_loc4(IsingLattice2D& lat) {
	//calculate the localization - absolute value of sz at each site, to the 4th power and averaged
	double final_result = 0.0;
	for (int i = 0; i < lat.get_Lx(); ++i) {
		double site_mag = 0.0;
		for (int j = 0; j < lat.get_Ly(); ++j) {
			site_mag += lat.get_spin(i, j);
		}
		site_mag /= lat.get_Ly();
		site_mag = site_mag*site_mag*site_mag*site_mag;
		final_result += site_mag;
	}
	return final_result / lat.get_Lx();
}

double GeneralLRW::calc_s1s2(IsingLattice2D& lat) {
	//calculate the spatial correlation function <S(i)S(i + 1)> (nearest neighbors only)
	double final_result = 0;
	double site_corr;
	int Lx = lat.get_Lx();
	for (int i = 0; i < Lx; ++i) {
		site_corr = 0;
		for (int j = 0; j < lat.get_Ly(); ++j) {
			site_corr += lat.get_spin(i, j)*lat.get_spin((i + 1)%Lx, j);
		}
		final_result += site_corr / lat.get_Ly();
	}
	return final_result / Lx;
}

double GeneralLRW::calc_sz_stagger(IsingLattice2D& lat) {
	//antiferromagnetic analogue of xmag.  At each imaginary time slice accumulate the sum over i of |(-1)^i * s_i(tau)|
	double final_result = 0.0;
	for (int j = 0; j < lat.get_Ly(); ++j) {
		double time_mag = 0.0;
		for (int i = 0; i < lat.get_Lx(); ++i) {
			time_mag += (2*(i/2) == i) ? lat.get_spin(i, j) : -lat.get_spin(i, j);
		}
		final_result += (time_mag > 0) ? time_mag / lat.get_Lx() : -time_mag / lat.get_Lx();
	}
	return final_result / lat.get_Ly();
}

double GeneralLRW::calc_action_slow(IsingLattice2D& lat) {
	//essentially, the time-averaged energy of a given set of quantum fluctuations
	double S = 0;
	int Lx = lat.get_Lx();
	int Ly = lat.get_Ly();
	int s1;
	for (int i = 0; i < Lx; ++i) {
		for (int j = 0; j < Ly; ++j) {
			s1 = lat.get_spin(i, j);
			for (int m = 0; m < Lx; ++m) {
				for (int n = 0; n < Ly; ++n) {
					S += s1 * lat.get_spin(m, n) * interactions[(m - i + Lx)%Lx][(n - j + Ly)%Ly];
				}
			}
		}
	}

	return 0.5*S;
}

double GeneralLRW::calc_point_action_slow(IsingLattice2D& lat, int x, int y){
	double S = 0;
	int Lx = lat.get_Lx();
	int Ly = lat.get_Ly();
	int s1 = lat.get_spin(x,y);
	for (int i = 0; i < Lx; ++i){
		for (int j = 0; j < Ly; ++j){
			S += s1 * lat.get_spin(i,j) * interactions[(i - x + Lx)%Lx][(j - y + Ly)%Ly];
		}
	}
	return S;
}

void GeneralLRW::step(IsingLattice2D& lat) {
	//1. determine a "seed" spin
	//2. go to all interacting spins and test to add to buffer and cluster
	//		- use long range procedure in Int. J. Mod. Phys. C 6, 359-370 (1995).
	//3. when buffer is empty, flip cluster
	cluster_size = 0;
	cluster_mag = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster[seed.x][seed.y] = 1;
	++cluster_size;
	cluster_mag += seed.proj;

	//2. add to buffer, then set loop to go through buffer
	test_spins(seed, lat);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins(seed, lat);
	}

	//3. flip cluster according to a probability check based on the magnetization change
	//prob is 1 if h*mag is negative, exp(-2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
	if (drand1_() < exp(-2 * h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster[x][y] == 1) {
					lat.flip_spin(x, y);
					cluster[x][y] = 0;
				}
			}
		}
	}
	mag -= (2.0*cluster_mag)/lat.get_N();

}

bool GeneralLRW::step_one_site(IsingLattice2D& lat, double * prev_action, thrust::host_vector<double>& corr_ref){
    //Build a cluster in only the time direction of one spatial site.  Then attempt to flip with Metropolis probability according to previous action
    //and newly calculated action.  Require GPU.  Mostly copied from step() function for full lattice clusters
    //return true if flipped, false if not
	cluster_size = 0;
	cluster_mag = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster[seed.x][seed.y] = 1;
	++cluster_size;
	cluster_mag += seed.proj;

	//2. add to buffer, then set loop to go through buffer
	timer.flag_start_time("Cluster Building");
	test_spins_one_site(seed, lat);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins_one_site(seed, lat);
	}
	timer.flag_end_time("Cluster Building");

	//3. flip cluster according to a probability check based on the magnetization change and interaction change
	//prob is 1 if( new_action - prev_action + 2*h*mag) is negative, exp(prev_action - new_action - 2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
    //
    //need to subtract off the action contribution from the seed site, because this portion of
    //detailed balance is already satisfied by the selection algorithm.
    double new_action = 0, new_t_cont = 0, prev_t_cont = 0;
	timer.flag_start_time("thrust variable creation");

	//single-site thrust vectors
    thrust::host_vector<double> thrustsite(lat.get_Ly());
    thrust::host_vector<double>& thrustsite_ref = thrustsite;
    thrust::host_vector<double> thrustint_site(lat.get_Ly());
    thrust::host_vector<double>& thrustint_site_ref = thrustint_site;

	//full lattice thrust vectors
	thrust::host_vector<double> thrustlat(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& thrustlat_ref = thrustlat;
	thrust::host_vector<double> thrustint(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& thrustint_ref = thrustint;
	lat.get_thrust_vector(thrustlat_ref);
	get_thrust_interactions(thrustint_ref);

	//throwaways for temporary correlation functions
	thrust::host_vector<double> site_corr_trash(lat.get_Ly());
	thrust::host_vector<double>& site_corr_trash_ref = site_corr_trash;
	thrust::host_vector<double> full_corr_trash(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& full_corr_trash_ref = full_corr_trash;

    for(int y = 0; y < lat.get_Ly(); ++y){
        thrustsite[y] = lat.get_spin(seed.x, y);
        thrustint_site[y] = interactions[0][y];
    }
	timer.flag_end_time("thrust variable creation");

	timer.flag_start_time("acceptance testing");
	prev_t_cont = cufft_calc_action(thrustsite_ref, thrustint_site_ref, site_corr_trash, 1, lat.get_Ly());

	//flip cluster to get proposed state
	for (int y = 0; y < lat.get_Ly(); ++y) {
		if (cluster[seed.x][y] == 1) {
			thrustlat[seed.x*lat.get_Ly() + y] *= -1;
			thrustsite[y] *= -1;
		}
	}
    new_action = cufft_calc_action(thrustlat_ref, thrustint_ref, full_corr_trash_ref, lat.get_Lx(), lat.get_Ly());
    new_t_cont = cufft_calc_action(thrustsite_ref, thrustint_site_ref, site_corr_trash, 1, lat.get_Ly());


                                
	if (drand1_() < exp((*prev_action - prev_t_cont) - (new_action - new_t_cont) - 2 * h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster[x][y] == 1) {
					lat.flip_spin(x, y);
					cluster[x][y] = 0;
				}
				//copy new correlation function
				corr_ref[x*lat.get_Ly() + y] = full_corr_trash[x*lat.get_Ly() + y];
			}
		}
		timer.flag_end_time("acceptance testing");
	    mag -= (2.0*cluster_mag)/lat.get_N();
        *prev_action = new_action;
        return true;
	}
    else{

		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
                cluster[x][y] = 0;
			}
		}
		timer.flag_end_time("acceptance testing");
        return false;
    }
}

bool GeneralLRW::step_one_site_prealloc(IsingLattice2D& lat, double * prev_action, thrust::host_vector<double>& corr_ref, 
		cufftHandle *full_forward_plan, cufftHandle *full_backward_plan, cufftHandle *onesite_forward_plan, cufftHandle *onesite_backward_plan,
		cudaError cuda_status, cufftDoubleReal *full_state_rs, cufftDoubleComplex *full_state_ft, cufftDoubleReal *onesite_state_rs, cufftDoubleComplex *onesite_state_ft){
    //Build a cluster in only the time direction of one spatial site.  Then attempt to flip with Metropolis probability according to previous action
    //and newly calculated action.  Require GPU.  Mostly copied from step() function for full lattice clusters
    //return true if flipped, false if not
	cluster_size = 0;
	cluster_mag = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster[seed.x][seed.y] = 1;
	++cluster_size;
	cluster_mag += seed.proj;

	//2. add to buffer, then set loop to go through buffer
	timer.flag_start_time("Cluster Building");
	test_spins_one_site(seed, lat);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins_one_site(seed, lat);
	}
	timer.flag_end_time("Cluster Building");

	//3. flip cluster according to a probability check based on the magnetization change and interaction change
	//prob is 1 if( new_action - prev_action + 2*h*mag) is negative, exp(prev_action - new_action - 2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
    //
    //need to subtract off the action contribution from the seed site, because this portion of
    //detailed balance is already satisfied by the selection algorithm.
    double new_action = 0, new_t_cont = 0, prev_t_cont = 0;
	timer.flag_start_time("thrust variable creation");

	//single-site thrust vectors
    thrust::host_vector<double> thrustsite(lat.get_Ly());
    thrust::host_vector<double>& thrustsite_ref = thrustsite;
    thrust::host_vector<double> thrustint_site(lat.get_Ly());
    thrust::host_vector<double>& thrustint_site_ref = thrustint_site;

	//full lattice thrust vectors
	thrust::host_vector<double> thrustlat(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& thrustlat_ref = thrustlat;
	thrust::host_vector<double> thrustint(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& thrustint_ref = thrustint;
	lat.get_thrust_vector(thrustlat_ref);
	get_thrust_interactions(thrustint_ref);

	//throwaways for temporary correlation functions
	thrust::host_vector<double> site_corr_trash(lat.get_Ly());
	thrust::host_vector<double>& site_corr_trash_ref = site_corr_trash;
	thrust::host_vector<double> full_corr_trash(lat.get_Lx()*lat.get_Ly());
	thrust::host_vector<double>& full_corr_trash_ref = full_corr_trash;

    for(int y = 0; y < lat.get_Ly(); ++y){
        thrustsite[y] = lat.get_spin(seed.x, y);
        thrustint_site[y] = interactions[0][y];
    }
	timer.flag_end_time("thrust variable creation");

	timer.flag_start_time("acceptance testing");
	prev_t_cont = cufft_calc_action_prealloc(thrustsite_ref, thrustint_site_ref, site_corr_trash, 1, lat.get_Ly(), 
						*onesite_forward_plan, *onesite_backward_plan, cuda_status, onesite_state_rs, onesite_state_ft);

	//flip cluster to get proposed state
	for (int y = 0; y < lat.get_Ly(); ++y) {
		if (cluster[seed.x][y] == 1) {
			thrustlat[seed.x*lat.get_Ly() + y] *= -1;
			thrustsite[y] *= -1;
		}
	}
    new_action = cufft_calc_action_prealloc(thrustlat_ref, thrustint_ref, full_corr_trash_ref, lat.get_Lx(), lat.get_Ly(), 
						*full_forward_plan, *full_backward_plan, cuda_status, full_state_rs, full_state_ft);
    new_t_cont = cufft_calc_action_prealloc(thrustsite_ref, thrustint_site_ref, site_corr_trash, 1, lat.get_Ly(),
						*onesite_forward_plan, *onesite_backward_plan, cuda_status, onesite_state_rs, onesite_state_ft);


                                
	if (drand1_() < exp((*prev_action - prev_t_cont) - (new_action - new_t_cont) - 2 * h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster[x][y] == 1) {
					lat.flip_spin(x, y);
					cluster[x][y] = 0;
				}
				//copy new correlation function
				corr_ref[x*lat.get_Ly() + y] = full_corr_trash[x*lat.get_Ly() + y];
			}
		}
		timer.flag_end_time("acceptance testing");
	    mag -= (2.0*cluster_mag)/lat.get_N();
        *prev_action = new_action;
        return true;
	}
    else{

		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
                cluster[x][y] = 0;
			}
		}
		timer.flag_end_time("acceptance testing");
        return false;
    }
}

bool GeneralLRW::step_one_site_prealloc_lowmem(IsingLattice2D& lat, double * prev_action, thrust::host_vector<double>& corr_ref, 
		thrust::host_vector<double>& thrustsite_ref,thrust::host_vector<double>& thrustint_site_ref,thrust::host_vector<double>& thrustlat_ref,
		thrust::host_vector<double>& thrustint_ref,thrust::host_vector<double>& site_corr_trash_ref,thrust::host_vector<double>& full_corr_trash_ref,
		cufftHandle *full_forward_plan, cufftHandle *full_backward_plan, cufftHandle *onesite_forward_plan, cufftHandle *onesite_backward_plan,
		cudaError cuda_status, cufftDoubleReal *full_state_rs, cufftDoubleComplex *full_state_ft, cufftDoubleReal *onesite_state_rs, cufftDoubleComplex *onesite_state_ft){
    //Build a cluster in only the time direction of one spatial site.  Then attempt to flip with Metropolis probability according to previous action
    //and newly calculated action.  Require GPU.  Mostly copied from step() function for full lattice clusters
	//use lower amount of memory by allocating temporary variables outside the loop and giving them new values in this function for each step
    //return true if flipped, false if not
	cluster_size = 0;
	cluster_mag = 0;
	//1. get seed
	spin seed;
    seed.x = (int) (lat.get_Lx() * drand1_());
    seed.y = (int) (lat.get_Ly() * drand1_());
    seed.proj = lat.get_spin(seed.x, seed.y);
    cluster[seed.x][seed.y] = 1;
	++cluster_size;
	cluster_mag += seed.proj;

	//2. add to buffer, then set loop to go through buffer
	timer.flag_start_time("Cluster Building");
	test_spins_one_site(seed, lat);
	while (buffer.size() != 0) {
		seed = buffer.back();
		cluster_mag += seed.proj;
		++cluster_size;
		buffer.pop_back();
		test_spins_one_site(seed, lat);
	}
	timer.flag_end_time("Cluster Building");

	//3. flip cluster according to a probability check based on the magnetization change and interaction change
	//prob is 1 if( new_action - prev_action + 2*h*mag) is negative, exp(prev_action - new_action - 2*h*mag) if h*mag is positive
	//satisfies detailed balance: this is the "acceptance probability", rest was "selection probability"
    //
    //need to subtract off the action contribution from the seed site, because this portion of
    //detailed balance is already satisfied by the selection algorithm.
    double new_action = 0, new_t_cont = 0, prev_t_cont = 0;
	lat.get_thrust_vector(thrustlat_ref);
	get_thrust_interactions(thrustint_ref);

    for(int y = 0; y < lat.get_Ly(); ++y){
        thrustsite_ref[y] = lat.get_spin(seed.x, y);
        thrustint_site_ref[y] = interactions[0][y];
    }

	timer.flag_start_time("acceptance testing");
	prev_t_cont = cufft_calc_action_prealloc(thrustsite_ref, thrustint_site_ref, site_corr_trash_ref, 1, lat.get_Ly(), 
						*onesite_forward_plan, *onesite_backward_plan, cuda_status, onesite_state_rs, onesite_state_ft);

	//flip cluster to get proposed state
	for (int y = 0; y < lat.get_Ly(); ++y) {
		if (cluster[seed.x][y] == 1) {
			thrustlat_ref[seed.x*lat.get_Ly() + y] *= -1;
			thrustsite_ref[y] *= -1;
		}
	}
    new_action = cufft_calc_action_prealloc(thrustlat_ref, thrustint_ref, full_corr_trash_ref, lat.get_Lx(), lat.get_Ly(), 
						*full_forward_plan, *full_backward_plan, cuda_status, full_state_rs, full_state_ft);
    new_t_cont = cufft_calc_action_prealloc(thrustsite_ref, thrustint_site_ref, site_corr_trash_ref, 1, lat.get_Ly(),
						*onesite_forward_plan, *onesite_backward_plan, cuda_status, onesite_state_rs, onesite_state_ft);


                                
	if (drand1_() < exp((*prev_action - prev_t_cont) - (new_action - new_t_cont) - 2 * h * cluster_mag)) {
		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
				if (cluster[x][y] == 1) {
					lat.flip_spin(x, y);
					cluster[x][y] = 0;
				}
				//copy new correlation function
				corr_ref[x*lat.get_Ly() + y] = full_corr_trash_ref[x*lat.get_Ly() + y];
			}
		}
		timer.flag_end_time("acceptance testing");
	    mag -= (2.0*cluster_mag)/lat.get_N();
        *prev_action = new_action;
        return true;
	}
    else{

		for (int x = 0; x < lat.get_Lx(); ++x) {
			for (int y = 0; y < lat.get_Ly(); ++y) {
                cluster[x][y] = 0;
			}
		}
		timer.flag_end_time("acceptance testing");
        return false;
    }
}


bool GeneralLRW::metropolis_step(IsingLattice2D& lat, double *action){
	int x = (int) (lat.get_Lx() * drand1_());
    int y = (int) (lat.get_Ly() * drand1_());
	double point_action = calc_point_action_slow(lat,x, y);
	if(drand1_() < exp(2.0*(point_action))){
		lat.flip_spin(x, y);
		*action -= 2*point_action;
		return true;
	}
	else{return false;}
}

void GeneralLRW::test_spins(spin seed, IsingLattice2D& lat){
	int i,j, k, l; //i and j is searched site distance, k and l are the corresponding lattice indices
	double randnum, totalint = 0; //random number, and running interaction sum
	std::vector<double>::iterator x_it_temp = site_sums.begin(), 
			x_it = site_sums.begin(), x_it_end = site_sums.end(), y_it, y_it_end; //x_it is 0 site separation, x_it_end is the end of the site vector
			//y_it will be the beginning of the interaction sum determined by the site search, and y_it_end will be the corresponding ending iterator
	spin add;
	//
	//	Outline of algorithm:
	//		1. get random number.  Find corresponding interaction sum by sum = -1/2 log(rand)
	//		2. search for totalint + randsum; first site_sums, then full interaction sum at the given site
	//		3. add corresponding spin if not in cluster and interaction with seed is favorable
	//		4. update totalint
	//
	//		5. end when totalint + randsum exceeds all interaction sums
	//
	while(x_it_temp != x_it_end){
		randnum = drand1_();
		totalint = totalint - 0.5*log(randnum);
		x_it_temp = std::lower_bound(x_it, x_it_end, totalint);
		if(x_it_temp != x_it_end){
			i = x_it_temp - x_it;
			y_it = interaction_sum[i].begin();
			y_it_end = interaction_sum[i].end();
			j = (std::lower_bound(y_it, y_it_end, totalint) - y_it);
			k = (seed.x + i) % lat.get_Lx();
			l = (seed.y + j) % lat.get_Ly();
			if (j < lat.get_Ly() && cluster[k][l] == 0 && interactions[i][j]*seed.proj*lat.get_spin(k, l) < 0){
				add.x = k;
				add.y = l;
				add.proj = lat.get_spin(k, l);
				buffer.push_back(add);
				cluster[k][l] = 1;
			}
		}
	}
}


void GeneralLRW::test_spins_one_site(spin seed, IsingLattice2D& lat){
    //same as test_spins(), but cluster is restricted to the original x position of the seed
	int i,j, k, l; //i and j is searched site distance, k and l are the corresponding lattice indices
	double randnum, totalint = 0; //random number, and running interaction sum
	std::vector<double>::iterator y_it, y_it_temp, y_it_end;//y_it will be the beginning of the interaction sum determined by the site search, and y_it_end will be the corresponding ending iterator
	spin add;
	//
	//	Outline of algorithm:
	//		1. get random number.  Find corresponding interaction sum by sum = -1/2 log(rand)
	//		2. search for totalint + randsum; first site_sums, then full interaction sum at the given site
	//		3. add corresponding spin if not in cluster and interaction with seed is favorable
	//		4. update totalint
	//
	//		5. end when totalint + randsum exceeds all interaction sums
	//

    i = seed.x;
    y_it = interaction_sum[0].begin();
	y_it_temp = y_it;
    y_it_end = interaction_sum[0].end();
    while(y_it_temp != y_it_end){
		randnum = drand1_();
		totalint = totalint - 0.5*log(randnum);
		y_it_temp = std::lower_bound(y_it, y_it_end, totalint);
		if(y_it_temp != y_it_end){
			j = (y_it_temp - y_it);
			l = (seed.y + j) % lat.get_Ly();
			if (j < lat.get_Ly() && cluster[i][l] == 0 && interactions[0][j]*seed.proj*lat.get_spin(i, l) < 0){
				add.x = i;
				add.y = l;
				add.proj = lat.get_spin(i, l);
				buffer.push_back(add);
				cluster[i][l] = 1;
			}
		}
    }
}
