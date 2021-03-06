#pragma once
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "IsingLattice2D.h"
extern "C" {
#include "random.h"
}

//classical monte carlo params for a number of different models and algorithms
//Models: n-vector models (ising, heisenberg), spinboson
//Lattices: chain, square, triangular, kagome
//Algorithms: wolff, generalizedwolff

const double PI = 3.141592653589793238462643;



struct spin_boson_params {
	double g, A0, delta, v, omega_c;
};

struct class_mc_params {
	int dim, rand_seed, eq_time, steps_per_measure, measures_per_dump, max_dumps;
	double beta, h, kT;
	std::vector<int> lengths;
	std::vector<double> spacings;
	std::vector<double> Js;
	std::string cutoff_type, lattice, alg;
	spin_boson_params sbparams;
	std::string to_string() {
		std::stringstream ss;
		ss << "Dimension: " << dim << "\n";
		ss << "Starting Lattice type: " << lattice << "\n";
		ss << "Dimension lengths: ";
		for (auto const& value : lengths) {
			ss << value << " ";
		}
		ss << "\n";
		ss << "Dimension spacings: ";
		for (auto const& value : spacings) {
			ss << value << " ";
		}
		ss << "\n";
		ss << "Cutoff: " << cutoff_type << "\n";
		ss << "J couplings: ";
		for (auto const& value : Js) {
			ss << value << " ";
		}
		ss << "\n";
		ss << "kT: " << kT << "\n";
		ss << "Beta: " << beta << "\n";
		ss << "H: " << h << "\n";
		ss << "Algorithm: " << alg << "\n";
		ss << "Random Seed: " << rand_seed << "\n";
		ss << "Equilibration time: " << eq_time << "\n";
		ss << "Steps Per Measurement: " << steps_per_measure << "\n";
		ss << "Measures Per Dump: " << measures_per_dump << "\n";
		ss << "Max Dumps: " << max_dumps << "\n";
		ss << "Spin Boson Params:\n";
		ss << "g: " << sbparams.g << "\n";
		ss << "A0: " << sbparams.A0 << "\n";
		ss << "Delta: " << sbparams.delta << "\n";
		ss << "V: " << sbparams.v << "\n";
		ss << "omega_c: " << sbparams.omega_c << "\n";

		return ss.str();
	}
};



struct class_mc_measurements {
	//vector for the names (steps, energies, mags, etc.)
	//parallel vector for the recorded values
	std::vector<std::string> names;
	std::vector<std::vector<double>> values;

	std::vector<std::string> func_names;
	std::vector<std::vector<std::vector<double>>> functions;

	std::vector<std::vector<int>> states;

	void write_results(std::ofstream * f, std::string type){
		//f must be open. type says "vals", "func_names[i]", or "states"
		if(f->is_open()){
			if(type.compare("vals")==0){
				for (int i = 0; i < names.size(); ++i) {
        			*f << names[i] ;
					for(int j = 0; j < values[i].size(); ++j){
						*f << "," << values[i][j];
					}
					*f << "\n";
    			}
			}
			else if(type.compare("states")==0){
				for(int i = 0; i < states.size(); ++i){
					*f << states[i][0];
					for(int j = 1; j < states[i].size(); ++j){
						*f << "," << states[i][j] ;
					}
					*f << "\n";
				}
			}
			else{
				bool found = false;
				for (int i = 0; i < func_names.size(); ++i) {
					if (type.compare(func_names[i]) == 0) {
						found = true;
						//functions[i] is a std::vector<std::vector<double>>
						//functions[i][j] is a std::vector<double> (jth measurement of function i)
						//functions[i][j][k] is a double
						for(int j = 0; j < functions[i].size(); ++j){
							*f << functions[i][j][0];
							for(int k = 1; k < functions[i][j].size(); ++k){
								*f << "," << functions[i][j][k];
							}
							*f << "\n";
						}
					}
				}
				if (!found) {
					std::cout << type << " is not a valid observable to write to file\n";
				}
			}
		}
		else{
			std::cout << "Error writing results: file not opened before calling function\n";
		}
	}

	std::vector<double> get_vals(std::string identifier) {
		int name_ind = -1;
		for (int i = 0; i < names.size(); ++i) {
			if (identifier.compare(names[i]) == 0) {
				name_ind = i;
			}
		}
		if (name_ind == -1) {
			std::cout << identifier << " is not a valid observable value\n";
			return {0};
		}
		else {
			return values[name_ind];
		}
	}

	void record(std::string name, double val) {
		bool found = false;
		for (int i = 0; i < names.size(); ++i) {
			if (name.compare(names[i]) == 0) {
				found = true;
				values[i].push_back(val);
			}
		}
		if (!found) {
			std::cout << name << " is not a valid observable\n";
		}
	}

	void record(std::string name, std::vector<double> function) {
		bool found = false;
		for (int i = 0; i < func_names.size(); ++i) {
			if (name.compare(func_names[i]) == 0) {
				found = true;
				//functions[i] is a std::vector<std::vector<double>>
				//functions[i][j] is a std::vector<double> (jth measurement of function i)
				//functions[i][j][k] is a double
				if (function.size() == functions[i][0].size()) {
					functions[i].push_back(function);
				}
                else if (functions[i][0].size() == 0){
                    for(int n = 0; n < function.size(); ++n){
                        functions[i][0].push_back(function[n]);
                    }
                }
				else { found = false; }
			}
		}
		if (!found) {
			std::cout << name << " is not a valid observable function, or there was a size mismatch between recordings\n";
		}
	}

	void record(IsingLattice2D& lat){
		states.push_back(lat.get_bool_spins());
	}


	void record(std::vector<double> vals) {
		if (vals.size() != names.size()) {
			std::cout << "Attempt to record failed: too many/few values recorded\n";
		}
		else {
			for (int i = 0; i < vals.size(); ++i) {
				values[i].push_back(vals[i]);
			}
		}
	}
};

void create_input();

bool str_is_equal(std::string, std::string);

std::string vec2str(std::vector<double> vec);
std::string vec2str(std::vector<int> vec);//i can't make template functions compile to save my life

void read_input_ising(std::ifstream*, class_mc_params*);

void read_input_spin_boson(std::ifstream*, spin_boson_params*);

void read_input_mpi(std::ifstream*, std::string*, std::vector<double>&);

void apply_spin_boson_params(class_mc_params*);

void write_outputs(int, std::vector<int>, std::vector<double>, std::vector<double>, std::vector<double>);

void write_outputs_var(std::ofstream*, class_mc_measurements);

void write_params(class_mc_params params);

void write_state(int, IsingLattice2D&);

void write_state_pbm(int, IsingLattice2D&);

void read_state_pbm(int, IsingLattice2D&);

bool isDirExist(const std::string& path);

bool makePath(const std::string& path);

double mean(std::vector<double> vals);

double error(std::vector<double> vals, double mean, int bins);

double bootstrap(std::vector<double>, int, std::string);
