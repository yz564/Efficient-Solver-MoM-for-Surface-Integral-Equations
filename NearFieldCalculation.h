//Calculating the near fields after solving the system, near singularity extracted
//Zhong Yang @ duke 2021

#pragma once
#include "meshes.h"


template <typename T>
void nearFieldHelper(Matrix<complex<T>, Dynamic, Dynamic>* I_solv, vector<T> observe, RWGbasis<T> * rwg, TriElement<T> * Tj, Matrix<complex<T>, Dynamic, Dynamic>* ans, const int row, bool IsSingular) {
	complex<T> eps = Tj->eps_out;
	complex<T> mu = Tj->mu_out;
	complex<T> k = (rwg->wf)*sqrt(eps*mu);
	complex<T> ik = k * complex<T>(0, -1);//-jk
	complex<T> Z = sqrt(mu / eps);
	complex<T> FactorF = complex<T>(1.0, 0) / (ik*ik);
	complex<T> jkZ = Z * k * complex<T>(0, 1);

	int M = Tj->iGPt.size();
	for (int m = 0; m < M; ++m) {
		vector<T> source = Tj->iGPt[m];
		complex<T> w1 = myReal2Complex(rwg->innerGaussWeight[m]);
		vector<complex<T> > R = myVecReal2Complex(myVecMinus(observe, source));
		complex<T> Rabs = myReal2Complex(myVecDistance(observe, source));
		vector<complex<T> > gradG;
		complex<T> g;
		if (IsSingular) {
			if (abs(Rabs) < 1e-15) {
				g = ik;
				gradG = myVecNumProd(R, complex<T>(0,1/3.0/PI<T>)*k*k*k);
			}
			else {
				g = (exp(ik * Rabs) - complex<T>(1.0, 0)) / Rabs;
				gradG = myVecNumProd(R, (complex<T>(1.0, 0) - (complex<T>(1.0, 0) - Rabs * ik)*exp(Rabs*ik)) / Rabs / Rabs / Rabs);
			}
		}
		else {
			g = exp(ik * Rabs) / Rabs;
			gradG = myVecNumProd(R, (-(complex<T>(1.0, 0) - Rabs * ik)*exp(Rabs*ik)) / Rabs / Rabs / Rabs);
		}
		for (int j = 0; j < 3; ++j) {
			int edgej = Tj->edgeId[j];
			complex<T> signj = myReal2Complex(Tj->edgeSign[j]);
			complex<T> lj = myReal2Complex(Tj->length[j]);
			complex<T> area_j = myReal2Complex(Tj->area);
			vector<complex<T> > rho_j = myVecReal2Complex(myVecMinus(source, Tj->nodes[j]));
			vector<complex<T> > pq = myVecReal2Complex(myVecMinus(observe, Tj->nodes[j]));
			complex<T> I = (*I_solv)(edgej);
			(*ans).block(row, 3, 1, 3) += (myVec2Eigen(myVecCross(myVecNumProd(rho_j,I), gradG))*w1*lj*signj / complex<T>(2.0, 0)); //H-field
			(*ans).block(row, 0, 1, 3) += (myVec2Eigen(myVecMinus(myVecNumProd(rho_j, g*I / complex<T>(2.0, 0)) , myVecNumProd(gradG,I*FactorF)))*w1*lj*signj*jkZ); //E-field
			if (IsSingular && (m == 0)) {
				vector<T> singular = Tj->graglia_1(j, observe);
				vector<T> V1_(3);
				V1_[0] = singular[0];
				V1_[1] = singular[1];
				V1_[2] = singular[2];
				vector<complex<T> > V1 = myVecReal2Complex(V1_);
				complex<T> I1 = myReal2Complex(singular[3]);
				vector<T> V2_(3);
				V2_[0] = singular[4];
				V2_[1] = singular[5];
				V2_[2] = singular[6];
				vector<complex<T> > V2 = myVecReal2Complex(V2_);
				(*ans).block(row, 3, 1, 3) += myVec2Eigen(myVecCross(myVecNumProd(pq,I), V2))*lj*signj/area_j/complex<T>(2.0, 0); //H-field
				(*ans).block(row, 0, 1, 3) += myVec2Eigen(myVecMinus(myVecNumProd(V1, I / complex<T>(2.0, 0)), myVecNumProd(V2, I*FactorF)))*lj*signj*jkZ / area_j; //E-field
			}
		}
	}
}


template <typename T>
Matrix<complex<T>, Dynamic, Dynamic>* EHfieldCal(RWGbasis<T> * rwg, Matrix<complex<T>, Dynamic, Dynamic>* I_solv, InputInfo<T> * param) {
	int num_elem = rwg->num_elem;
	int N_ = param->receiver_numbers;
	vector<int> N(3, 1);
	vector<T> step(3, 0);
	if (param->receiver_start[0] != param->receiver_end[0]) {
		N[0]=N_;
		step[0] = (param->receiver_end[0] - param->receiver_start[0]) / (T)(N_ - 1);
	}
	if (param->receiver_start[1] != param->receiver_end[1]) {
		N[1] = N_;
		step[1] = (param->receiver_end[1] - param->receiver_start[1]) / (T)(N_ - 1);
	}
	if (param->receiver_start[2] != param->receiver_end[2]) {
		N[2] = N_;
		step[2] = (param->receiver_end[2] - param->receiver_start[2]) / (T)(N_ - 1);
	}
	Matrix<complex<T>, Dynamic, Dynamic> * ans = new Matrix<complex<T>, Dynamic, Dynamic>(N[0] * N[1] * N[2], 6);
	(*ans).setZero();
	int count = 0;
	std::ofstream outfile;
	outfile.open("receivers.txt", std::ios_base::trunc);
	for (int ix = 0; ix < N[0]; ++ix) {
		vector<T> receiver;
		receiver.push_back(param->receiver_start[0]+ix*step[0]);
		for (int iy = 0; iy < N[1]; ++iy) {
			receiver.push_back(param->receiver_start[1] + iy * step[1]);
			for (int iz = 0; iz < N[2]; ++iz) {
				receiver.push_back(param->receiver_start[2] + iz * step[2]);
				
				outfile<<"receiver: "<< receiver[0] << ", " << receiver[1] << ", " << receiver[2] << endl;
				int j;
				Matrix<complex<T>, Dynamic, Dynamic> tmp(num_elem, 6);
				tmp.setZero();
				#pragma omp parallel for private(j)
				for (j = 0; j<num_elem; ++j) {
					TriElement<T> * Tj = rwg->triangles[j];
					T distance = myVecDistance(receiver, Tj->iGPt[0]);
					if (distance < rwg->singularThreshold * Tj->length[0]) {
						nearFieldHelper<T>(I_solv, receiver, rwg, Tj, &tmp, j, true);
					}
					else {
						nearFieldHelper<T>(I_solv, receiver, rwg, Tj, &tmp, j, false);
					}
				}
				(*ans).row(count) = tmp.colwise().sum();
				count++;
				receiver.pop_back();
			}
			receiver.pop_back();
		}
		receiver.pop_back();
	}
	*ans = (*ans) / myReal2Complex(-4.0 * PI<T>);
	outfile.close();
	return ans;
}