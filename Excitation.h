//Calculating the excitation vector (the righ hand side), near singularity extracted.
//Only the M dipole calculation is implemented
//Zhong Yang @ duke 2021

#pragma once
#include "meshes.h"

template <typename T>
void rhsHelper_K_operator(vector<complex<T>> A, vector<T> source, RWGbasis<T> * rwg, TriElement<T> * Ti, Matrix<complex<T>, Dynamic, Dynamic>* ans, bool IsSingular) {
	complex<T> eps = Ti->eps_out;
	complex<T> mu = Ti->mu_out;
	complex<T> k = (rwg->wf)*sqrt(eps*mu);
	complex<T> ik = k * complex<T>(0, -1.0);//-jk
	//complex<T> Z(sqrt(mu / eps), 0);
	//complex<T> FactorF = 1;
	//FactorF /= (ik*ik);

	int M = Ti->oGPt.size();
	for (int m = 0; m < M; ++m) {
		vector<T> observe = Ti->oGPt[m];
		complex<T> w1 = myReal2Complex(rwg->outerGaussWeight[m]);
		vector<complex<T> > R = myVecReal2Complex(myVecMinus(observe, source));
		//complex<T> Rabs = myVecNorm(R); 
		complex<T> Rabs = myReal2Complex(myVecDistance(observe, source));
		vector<complex<T> > gradG;
		if (IsSingular) {
			gradG = myVecNumProd(R,(complex<T>(1.0, 0) -(complex<T>(1.0,0) -  Rabs*ik)*exp(Rabs*ik))/ Rabs/Rabs/Rabs);
		}
		else {
			gradG = myVecNumProd(R, (- (complex<T>(1.0, 0) - Rabs*ik )*exp(Rabs*ik)) / Rabs / Rabs / Rabs);
		}
		for (int i = 0; i < 3; ++i) {
			int edgei = Ti->edgeId[i];
			complex<T> signi = myReal2Complex(Ti->edgeSign[i]);
			complex<T> li = myReal2Complex(Ti->length[i]);
			complex<T> area_i= myReal2Complex(Ti->area);
			vector<complex<T> > rho_i= myVecReal2Complex(myVecMinus(observe, Ti->nodes[i]));
			vector<complex<T> > pq = myVecReal2Complex(myVecMinus(source, Ti->nodes[i]));
			(*ans)(edgei) += (myVecDot(rho_i,myVecCross(A,gradG))*w1*li*signi / complex<T>(2.0,0));
			if (IsSingular && (m == 0)) {
				vector<T> singular = Ti->graglia_1(i, source);
				vector<T> V2_(3);
				V2_[0] = singular[4];
				V2_[1] = singular[5];
				V2_[2] = singular[6];
				vector<complex<T> > V2 = myVecReal2Complex(V2_);
				(*ans)(edgei) += (li*signi*myVecDot(A, myVecCross(pq,V2)) / area_i / complex<T>(2.0, 0));
			}
		}
	}
}

//template <typename T>
//class EMSource {
//public:
//	T magnitude;
//	T phase;
//	T theta;
//	T phi;
//	virtual Matrix<complex<T>, Dynamic, Dynamic> * testEFieldCal(RWGbasis<T> * rwg)=0;
//	virtual Matrix<complex<T>, Dynamic, Dynamic> * testHFieldCal(RWGbasis<T> * rwg) = 0;
//
//};

template <typename T>
//class Mdipole: public EMSource<T> {
class Mdipole {
public:
	T magnitude;
	T phase;
	T theta;
	T phi;
	vector<T> location;
	Mdipole(T A, T p, T t_, T p_, vector<T> loc) : magnitude(A), phase(p), theta(t_), phi(p_), location(loc) {}
	Matrix<complex<T>, Dynamic, Dynamic> * testEFieldCal(RWGbasis<T> * rwg) {
		complex<T> jphase(0, phase);
		complex<T> M_ = magnitude * exp(jphase);
		vector<complex<T> > M(3);
		M[0] = M_*sin(theta)*cos(phi);
		M[1] = M_ * sin(theta)*sin(phi);
		M[2] = M_ * cos(theta);
		int num_edge = rwg->num_edge;
		int num_elem = rwg->num_elem;
		Matrix<complex<T>, Dynamic, Dynamic>* ans = new Matrix<complex<T>, Dynamic, Dynamic>(num_edge, 1);
		(*ans).setZero();
		int i=0;
		#pragma omp parallel for private(i)
		for (i = 0; i < num_elem; ++i) {
			TriElement<T> * Ti = rwg->triangles[i];
			T distance = myVecDistance(location, Ti->oGPt[0]);
			if (distance < rwg->singularThreshold * Ti->length[0]) {
				rhsHelper_K_operator<T>(M, location, rwg, Ti, ans, true);
			}
			else {
				rhsHelper_K_operator<T>(M, location, rwg, Ti, ans, false);
			}
		}
		*ans = (*ans) / myReal2Complex(4.0 * PI<T>);
		return ans;
	}
};