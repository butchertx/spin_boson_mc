#include "class_mc_io.h"
#include <cmath>
#include <iostream>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
	struct _stat info;
	if (_stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & _S_IFDIR) != 0;
#else 
	struct stat info;
	if (stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& path)
{
#if defined(_WIN32)
	int ret = _mkdir(path.c_str());
#else
	mode_t mode = 0755;
	int ret = mkdir(path.c_str(), mode);
#endif
	if (ret == 0)
		return true;

	switch (errno)
	{
	case ENOENT:
		// parent didn't exist, try to create it
	{
		int pos = path.find_last_of('/');
		if (pos == std::string::npos)
#if defined(_WIN32)
			pos = path.find_last_of('\\');
		if (pos == std::string::npos)
#endif
			return false;
		if (!makePath(path.substr(0, pos)))
			return false;
	}
	// now, try to create again
#if defined(_WIN32)
	return 0 == _mkdir(path.c_str());
#else 
	return 0 == mkdir(path.c_str(), mode);
#endif

	case EEXIST:
		// done!
		return isDirExist(path);

	default:
		return false;
	}
}

void create_input() {

}

bool str_is_equal(std::string str1, std::string str2) {
	return str1.compare(str2) == 0;
}

/*Need to read many params
list order here
*/
void read_input_ising(std::ifstream* file_p, class_mc_params* params) {
	std::string line;
	int dummyi;
	double dummyd;
	std::stringstream iss;

	if (file_p->is_open()) {
		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#lattice parameters

		std::getline(*file_p, line);
		iss << line;//#dimension
		iss >> params->dim;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#type
		iss >> params->lattice;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#lengths
		for (int i = 0; i < params->dim; ++i) {
			iss >> dummyi;
			params->lengths.push_back(dummyi);
		}
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#spacings
		for (int i = 0; i < params->dim; ++i) {
			iss >> dummyd;
			params->spacings.push_back(dummyd);
		}
		iss.str("");

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#blank

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#Model Parameters

		std::getline(*file_p, line);
		iss << line;//#cutoff type
		iss >> params->cutoff_type;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Jcouples
		for (int i = 0; i < params->dim; ++i) {
			iss >> dummyd;
			params->Js.push_back(dummyd);
		}
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#beta
		iss >> params->beta;
		iss.str("");

		params->kT = 1 / params->beta;

		std::getline(*file_p, line);
		iss << line;//#h
		iss >> params->h;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#blank

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#Algorithm Parameters

		std::getline(*file_p, line);
		iss << line;//#algorithm
		iss >> params->alg;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#random seed
		iss >> params->rand_seed;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#blank

		std::getline(*file_p, line);
		iss << line;
		iss.str("");//#Simulation Parameters

		std::getline(*file_p, line);
		iss << line;//#Metropolis steps
		iss >> params->metropolis_steps;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Onesite Steps
		iss >> params->onesite_steps;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Wolff Steps
		iss >> params->wolff_steps;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Kept Measures
		iss >> params->kept_measures;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Throwaway Measures
		iss >> params->throwaway_measures;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#Parallel Tempering Steps
		iss >> params->ptemp_steps;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#blank
		iss.str("");
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
}

void read_input_spin_boson(std::ifstream* file_p, spin_boson_params* params) {
	std::string line;
	int dummyi;
	double dummyd;
	std::stringstream iss;

	if (file_p->is_open()) {
		std::getline(*file_p, line);
		iss << line;//#Spin Boson Params
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#g
		iss >> params->g;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#A0
		iss >> params->A0;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#delta
		iss >> params->delta;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#v
		iss >> params->v;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#omega_c
		iss >> params->omega_c;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#blank
		iss.str("");
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
}

void read_input_mpi(std::ifstream* file_p, std::string* param_name, std::vector<double>& vals) {
	std::string line;
	int dummyi;
	double dummyd;
	std::stringstream iss;

	if (file_p->is_open()) {
		std::getline(*file_p, line);
		iss << line;//#MPI params
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#param name
		iss >> *param_name;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#num_vals
		iss >> dummyi;
		iss.str("");

		std::getline(*file_p, line);
		iss << line;//#vals
		for (int i = 0; i < dummyi; ++i){
            iss >> dummyd;
			vals.push_back(dummyd);
		}
		iss.str("");
	}
	else {
		std::cout << "Error: input file not opened\n";
	}
}


void apply_spin_boson_params(class_mc_params* params) {
	double tc = params->beta / ((double)params->lengths[1]);
	params->spacings[1] = tc;
	params->Js[0] = tc*params->Js[0];
	params->Js[1] = -0.5 * log(tanh(tc * params->sbparams.delta));
	params->h = tc*params->h;
}

std::string vec2str(std::vector<double> vec) {
	std::stringstream ss;
	for (int i = 0; i < vec.size() - 1; ++i) {
		ss << vec[i] << ", ";
	}
	ss << vec.back();
	return ss.str();
}

std::string vec2str(std::vector<int> vec) {
	std::stringstream ss;
	for (int i = 0; i < vec.size() - 1; ++i) {
		ss << vec[i] << ", ";
	}
	ss << vec.back();
	return ss.str();
}

void write_outputs(int dump_num, std::vector<int> steps, std::vector<double> times, std::vector<double> record1, std::vector<double> record2) {
	//char* dump_path = "./dump";
	makePath("./dump");
	char dump_name[100];
	sprintf(dump_name, "dump/dump%d.csv", dump_num);
	std::ofstream file;
	file.open(dump_name);
	file << "Steps," << vec2str(steps) << "\n";
	file << "Times," << vec2str(times) << "\n";
	file << "<Sz>," << vec2str(record1) << "\n";
	file << "<Sx>," << vec2str(record2) << "\n";
	//file << "corr_t," << vec2str(record3) << "\n";
	file.close();
}

void write_outputs_var(std::ofstream * file_p, class_mc_measurements results) {
	if(file_p->is_open()){
		for (std::string name : results.names) {
			*file_p << name << "," << vec2str(results.get_vals(name)) << "\n";
		}
	}
}

void write_state(int state_num, IsingLattice2D& lat) {
	//char* dump_path = "./dump";
	//first line is Lx, Ly
	//following lines are 1's and 0's for up and down, with spaces in between
	makePath("./dump");
	char dump_name[100];
	sprintf(dump_name, "dump/state%d.csv", state_num);
	std::ofstream file;
	file.open(dump_name);
	file << lat.get_Lx() << " " << lat.get_Ly() << "\n";
	file << lat.to_string();
	file.close();

}

void write_state_pbm(int state_num, IsingLattice2D& lat){
	makePath("./dump");
	char dump_name[100];
	sprintf(dump_name, "dump/state%04d.pbm", state_num);
	std::ofstream file;
	file.open(dump_name);
	file << "P1 " << lat.get_Lx() << " " << lat.get_Ly() << "\n";
	file << lat.to_bitmap();
	file.close();
}

void read_state_pbm(int state_num, IsingLattice2D& lat){
	std::ifstream file;
	char filename[100];
	int Lx,Ly, dummyi;
	std::stringstream iss;
	std::string line;
	sprintf(filename, "dump/state%d.pbm", state_num);
	file.open(filename);

	std::getline(file, line);
	iss << line;//Lx and Ly
	iss >> Lx;
	iss >> Ly;
	iss.str("");
	if(Lx == lat.get_Lx() && Ly == lat.get_Ly()){
		for(int y = 0; y < Ly; ++y){
			std::getline(file, line);
			iss << line;
			for(int x = 0; x < Lx; ++x){
				iss >> dummyi;
				lat.set_spin(x,y,dummyi);
			}
			iss.str("");
		}
	}
	else{
		std::cout << "Warning: read state dimension mismatch.  Starting with default lattice\n";
	}

	file.close();
}

//find the mean of a given list of values
double mean(std::vector<double> vals) {
	double sum = 0.0;
	for (int i = 0; i < vals.size(); ++i) {
		sum += vals[i];
	}
	return sum / vals.size();
}

//find the error via the bin technique using a specified number of bins
double error(std::vector<double> vals, double mean, int bins) {
	if (bins < 2) {
		return 0.0;
	}
	else {
		int NMC = bins;//Nb is bin size, NMC is number of bins
		int Nb = vals.size() / NMC;
		std::vector<double> avgs(NMC);
		for (int i = 0; i < NMC; ++i) {
			double avg = 0;
			for (int j = 0; j < Nb; ++j) {
				avg += vals[Nb*i + j];
			}
			avg = avg / Nb;
			avgs[i] = avg;
		}
		double std_dev = 0;
		for (int i = 0; i < avgs.size(); ++i) {
			std_dev += (avgs[i] - mean)*(avgs[i] - mean);
		}
		std_dev = sqrt(std_dev / NMC / (NMC - 1));
		return std_dev;
	}
}

double jackknife_binder(std::vector<double> mag2, std::vector<double> mag4, int bins) {
	double binder_avg = 0;
	for (int i = 0; i < mag2.size(); ++i) {
		binder_avg += mag4[i] / mag2[i] / mag2[i];
	}
	binder_avg /= mag2.size();
	std::vector<double> b_jacks(bins, 0.0);
	for (int b = 0; b < bins; ++b) {
		for (int i = 0; i < mag2.size() / bins; ++i) {

		}
	}
	return 0;
}

double bootstrap(std::vector<double> x, int n_boot, std::string func) {
	//compute the error in <x^2> - <x>^2 using the bootstrap method
	//use n_boot as the number of bootstrap bins and a random seed
	std::vector<double> boots(n_boot);
	rand_init_(&n_boot);
	double xB, xB_sq, x_temp;
	int N = x.size();
	for (int boot = 0; boot < n_boot; ++boot) {
		//find a bootstrap value for x^B_boot and x^2 ^B_boot
		xB = 0;
		xB_sq = 0;
		for (int i = 0; i < N; ++i) {
			x_temp = x[(int)(drand1_()*N)];
			xB += x_temp;
			xB_sq += x_temp*x_temp;
		}
		xB /= N;
		xB_sq /= N;
		//find the values for f^B_boot
		if (func.compare("binder") == 0) {
			//for this case, input x should be m^2 measurements
			boots[boot] = xB_sq / (xB*xB);
		}
		else if (func.compare("susceptibility") == 0) {
			boots[boot] = xB_sq - xB*xB;
		}
		else {
			std::cout << "Error: invalid bootstrap function option \"" << func << "\". Returning 0 for error estimate.\n";
			return 0;
		}
	}
	//average the values for f^B_boot to get susceptibility
	double av = mean(boots);
	//std::cout << "Bootstrap mean: " << av << "\n";

	//std dev for f^B is std dev for f
	double boots_sq = 0;
	for (int i = 0; i < n_boot; ++i) {
		boots_sq += boots[i] * boots[i];
	}
	boots_sq /= n_boot;
	return sqrt(boots_sq - av*av);
}
