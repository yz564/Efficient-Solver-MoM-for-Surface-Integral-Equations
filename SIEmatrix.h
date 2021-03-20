//Calculating the impedance matrix for PEC (EFIE), near singularity extracted
//Zhong Yang @ duke 2021
#ifndef SIEMATRIX_H// include guard
#define SIEMATRIX_H
#include "meshes.h"
#include <omp.h>


template <typename T>
void ImpedanceMatrixHelper_L_operator(Matrix<complex<T>, Dynamic, Dynamic> *ans, RWGbasis<T>* rwg, TriElement<T> * Ti, TriElement<T> * Tj, bool IsSingular) {
	complex<T> eps = Ti->eps_out;
	complex<T> mu = Ti->mu_out;
	complex<T> k = myReal2Complex(rwg->wf)*sqrt(eps*mu);
	complex<T> ik = k * complex<T>(0, -1.0);//-jk
	complex<T> Z=sqrt(mu / eps);
	complex<T> FactorF = complex<T>(1.0, 0) / (ik*ik);
	complex<T> jkZ = Z* k * complex<T>(0, 1);

	int M = Ti->oGPt.size();
	int N = Tj->iGPt.size();
	for (int m = 0; m < M; ++m) {
		vector<T> observe = Ti->oGPt[m];
		complex<T> w1 = myReal2Complex(rwg->outerGaussWeight[m]);
		for (int n = 0; n < N; ++n) {
			vector<T> source = Tj->iGPt[n];
			complex<T> w2 = myReal2Complex(rwg->innerGaussWeight[n]);
			complex<T> w1w2 = w1 * w2;
			complex<T> R = myReal2Complex(myVecDistance(observe, source));
			complex<T> g;
			if (IsSingular) {
				if (abs(R)<1e-15) {
					g = ik;
				}
				else {
					g = (exp(ik * R) - complex<T>(1.0, 0)) / R;
				}
			}
			else {
				g = exp(ik * R) / R;
			}
			for (int i = 0; i < 3; ++i) {
				int edgei = Ti->edgeId[i];
				complex<T> signi = myReal2Complex(Ti->edgeSign[i]);
				complex<T> li = myReal2Complex(Ti->length[i]);
				vector<T> rho_i = myVecMinus(observe, Ti->nodes[i]);
				for (int j = 0; j < 3; ++j) {
					int edgej = Tj->edgeId[j];
					complex<T> signj = myReal2Complex(Tj->edgeSign[j]);
					complex<T> sign = signi * signj;
					complex<T> lj = myReal2Complex(Tj->length[j]);
					complex<T> lilj = li * lj;
					vector<T> rho_j = myVecMinus(source, Tj->nodes[j]);
					complex<T> area_j = myReal2Complex(Tj->area);
					complex<T> FactorA = myReal2Complex(myVecDot(rho_i, rho_j) / 4.0);
					(*ans)(edgei, edgej) += (w1w2*lilj*sign*(FactorA + FactorF)*g)*jkZ;
					if (IsSingular && (n == 0)) {
						vector<T> singular = Tj->graglia_1(j, observe);
						vector<T> V1(3);
						V1[0] = singular[0];
						V1[1] = singular[1];
						V1[2] = singular[2];
						complex<T> I1 = myReal2Complex(singular[3]);
						(*ans)(edgei, edgej) += (w1*lilj*sign*(I1*FactorF+ myReal2Complex(myVecDot(V1, rho_i))/ complex<T>(4.0, 0))/ area_j)*jkZ;
					}
				}
			}
		}
	}
}

template <typename T>
Matrix<complex<T>, Dynamic, Dynamic> * SIEmatrixCal_PEC(InputInfo<T> * param, RWGbasis<T>* rwg) {
	rwg->GaussPoints(param->innerGaussDegree, param->outerGaussDegree);//inner (source) Gauss integral degree and outer (testing) Gauss integral degree
	int num_edge = rwg->num_edge;
	int num_elem = rwg->num_elem;


	Matrix<complex<T>, Dynamic, Dynamic> * ans = new Matrix<complex<T>, Dynamic, Dynamic>(num_edge, num_edge);
	(*ans).setZero();
	
	int i, j;
#pragma omp parallel for private(i,j)
	for (i = 0; i < num_elem; ++i) {
		TriElement<T> * Ti = rwg->triangles[i];
		for (j = 0; j < num_elem; ++j) {
			TriElement<T> * Tj = rwg->triangles[j];
			T distance = myVecDistance(Ti->oGPt[0], Tj->iGPt[0]);
			if (distance < rwg->singularThreshold* Ti->length[0]) {
				ImpedanceMatrixHelper_L_operator<T>(ans, rwg, Ti, Tj, true);
			}
			else {
				ImpedanceMatrixHelper_L_operator<T>(ans, rwg, Ti, Tj, false);
			}
		}
	}
	*ans = (*ans) / myReal2Complex(4.0 * PI<T>);
	return ans;
}




#endif