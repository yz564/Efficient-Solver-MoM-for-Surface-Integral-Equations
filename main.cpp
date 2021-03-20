//This project is an efficient EFIE solver for PEC objects
//Calculating Z matrix is just a little bit slower than FEKO at the same CPU usage. Set more threads can speed up!
//Near singularity extraction included, the accuracy for near field problem is ensured
//Zhong Yang @ duke 2021
#include "parameters.h"
#include "MeshProcessing.h"
#include "SIEmatrix.h"
#include "Excitation.h"
#include "NearFieldCalculation.h"
#include <stdlib.h>
#include <ctime>
#include <omp.h>
#include <ppl.h>
#include <array>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include "write_bmm_ver01.h"
using namespace std;
using namespace concurrency;


template<typename Derived>
inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
{
	return ((x - x).array() == (x - x).array()).all();
}

template<typename Derived>
inline bool is_nan(const Eigen::MatrixBase<Derived>& x)
{
	return ((x.array() == x.array())).all();
}

void timer_helper(const string& message, ofstream & logfile, vector<clock_t>& timer) {
	cout << message << (clock() - timer.back()) / (double)CLOCKS_PER_SEC << " seconds" << endl;
	logfile << message << (clock() - timer.back()) / (double)CLOCKS_PER_SEC << " seconds" << endl;
	timer.push_back(clock());
};

template <typename T>
void pureSIE(){
	ofstream logfile;
	logfile.open("logfile.txt", std::ios_base::trunc);
	vector<clock_t> timer;
	cout<<"==================start=================="<<endl;
	logfile << "==================start==================" << endl;
	timer.push_back(clock());
	string ConfigFile = "Config.txt";
	
	InputInfo<T> * param = ReadInput<T>(ConfigFile);

	omp_set_num_threads(param->allowedthreads);
	Eigen::setNbThreads(param->allowedthreads);
	string meshfile = param->mesh_file;
	Mesh<T> * SIEmesh = LoadMesh<T>(meshfile);
	
	timer_helper("Load the mesh takes ", logfile, timer);

	RWGbasis<T>*rwg = PrecalRWGbasis<T>(SIEmesh,param);
	delete SIEmesh;

	timer_helper("Processing mesh takes ", logfile, timer);
	int num_edge = rwg->num_edge;
	cout << "the number of elements is " << rwg->num_elem << endl;
	cout << "the number of edges is " << num_edge << endl;
	cout << "the average edge length is " << rwg->AvgEdgeLength << endl;
	logfile << "the number of elements is " << rwg->num_elem << endl;
	logfile << "the number of edges is " << num_edge << endl;
	logfile << "the average edge length is " << rwg->AvgEdgeLength << endl;

	cout << "Calculating the Z matrix, please be patient :)" << endl;

	Matrix<complex<T>, Dynamic, Dynamic> * IMPMAT = SIEmatrixCal_PEC<T>(param,rwg);

	cout << "Exporting the impedance matrix Z :)" << endl;
	string ff = "Zmatrix_" + meshfile.substr(0, meshfile.length() - 4) + ".bmm";
	const char * fpt = ff.c_str();
	write_bmm(*IMPMAT, fpt);


	timer_helper("Impedance matrix calculation takes ", logfile, timer);

	vector<T> magnitude=param->source_magnitude;
	vector<T> phase=param->source_phase;
	vector<T> theta=param->source_theta;
	vector<T> phi=param->source_phi;
	Matrix<complex<T>, Dynamic, Dynamic> * RHS = new Matrix<complex<T>, Dynamic, Dynamic>(num_edge, 1);
	(*RHS).setZero();
	for (int i = 0; i < param->source_num; ++i) {
		if (param->source_type[i] == 1) {
			Mdipole<T> * source = new Mdipole<T>(magnitude[i], phase[i], theta[i], phi[i], param->source_location[i]);
			(*RHS) += *(source->testEFieldCal(rwg));
			delete source;
		}
	}

	timer_helper("Excitation vector calculation takes ", logfile, timer);
	cout << "Exporting the right hand side vector :)" << endl;
	string ff1 = "RHS_" + meshfile.substr(0, meshfile.length() - 4) + ".bmm";
	const char * fpt1 = ff1.c_str();
	write_bmm(*RHS, fpt1);


	cout << "Solving the system using the Eigen LU solver :)" << endl;
	Matrix<complex<T>, Dynamic, Dynamic> I_solv(num_edge, 1);
	I_solv = (*IMPMAT).lu().solve(*RHS);


	timer_helper("Solving the fully dense system takes ", logfile, timer);

	cout << "Calculating the near fields Escat and Hscat :)" << endl;
	Matrix<complex<T>, Dynamic, Dynamic>* nearfield = EHfieldCal(rwg, &I_solv, param);

	timer_helper("Calculating E and H near fields takes ", logfile, timer);

	cout << "Exporting the near fields Escat and Hscat :)" << endl;
	string ff2 = "ScatterFields_" + meshfile.substr(0, meshfile.length() - 4) + ".bmm";
	const char * fpt2 =ff2.c_str();
	write_bmm(*nearfield, fpt2);

	cout << "Total time: " << (clock() - timer.front()) / (double)CLOCKS_PER_SEC << " seconds" << endl;
	logfile << "Total time: " << (clock() - timer.front()) / (double)CLOCKS_PER_SEC << " seconds" << endl;
	logfile.close();
}

int main() {
	pureSIE<double>();
	return 0;
}