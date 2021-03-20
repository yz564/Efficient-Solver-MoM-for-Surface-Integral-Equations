//All in elements. the class TriElement includes all infomation for the calculation, its method graglia_1 is for the near singularity extraction
//Zhong Yang @ duke 2021
#ifndef MESHES_H// include guard
#define MESHES_H
#include <iostream>
#include <fstream>
#include<string>
#include <cmath>
#include "GaussPoints.h"
#include <assert.h>
#include "parameters.h"
using namespace std;
using namespace Eigen;

template<typename T>
class TriElement {
public:
	vector<vector<T> > nodes;
	T area;
	vector<T> length;
	vector<int> edgeId;
	vector<T> edgeSign;
	complex<T> eps_out;
	complex<T> mu_out;
	complex<T> eps_in;
	complex<T> mu_in;
	T u3, v3;
	vector<T> u_hat;
	vector<T> v_hat;
	vector<T> w_hat;
	vector<vector<T> > m_hat;
	vector<vector<T> > iGPt;
	vector<vector<T> > oGPt;


	TriElement(const vector<T>& n1, const vector<T>& n2, const vector<T>& n3);
	void initialMaterials(T wf, T eps_r1, T mu_r1, T cond1, T losstangent1, T eps_r2, T mu_r2, T cond2, T losstangent2);
	void initialGaussPoints(int indegree, int outdegree);
	vector<T> graglia_1(int freevertex, const vector<T>& observePoint);
};

template<typename T>
class Mesh {
public:
	int num_elem;
	int num_node;
	int elem_type;
	vector<vector<int>> elems;
	vector<int> attribs;
	vector<vector<T>> nodes;

	Mesh(int Ne, int Nn, int Et) :num_elem(Ne), num_node(Nn), elem_type(Et){}


};

template<typename T>
class RWGbasis{
public:
	int num_edge;
	int num_elem;
	T singularThreshold;
	T wf;//2*PI*freq
	T AvgEdgeLength;
	vector<TriElement<T>*> triangles;
	vector<T> innerGaussWeight;
	vector<T> outerGaussWeight;


	
	~RWGbasis() {
		for (int i = 0; i < triangles.size(); ++i) {
			delete triangles[i];
		}
	}
	void findSharedEdges(const Mesh<T> *, const InputInfo<T> *);
	void GaussPoints(int indegree, int outdegree);

};



//member method implementation:
template<typename T>
TriElement<T>::TriElement(const vector<T>& n1, const vector<T>& n2, const vector<T>& n3) {
	nodes.resize(3);
	nodes[0] = n1;
	nodes[1] = n2;
	nodes[2] = n3;
	edgeId.resize(3);
	edgeSign.resize(3);
	vector<T> p21 = myVecMinus<T>(n2, n1);
	vector<T> p31 = myVecMinus<T>(n3, n1);
	vector<T> n_ = myVecCross<T>(p21, p31);
	T n_abs = myVecNorm<T>(n_);
	w_hat = myVecNumDiv<T>(n_, n_abs);
	area = n_abs / 2;
	vector<T> p32 = myVecMinus<T>(n3, n2);
	length.resize(3);
	length[0] = myVecNorm<T>(p32);
	length[1] = myVecNorm<T>(p31); //p13
	length[2] = myVecNorm<T>(p21);
	m_hat.resize(3);
	m_hat[0] = myVecCross<T>(myVecNumDiv<T>(p32, length[0]), w_hat);
	m_hat[1] = myVecCross<T>(w_hat, myVecNumDiv<T>(p31, length[1]));
	m_hat[2] = myVecCross<T>(myVecNumDiv<T>(p21, length[2]), w_hat);
	u_hat = myVecNumDiv<T>(p21, length[2]);
	v_hat = myVecCross<T>(w_hat, u_hat);
	u3 = myVecDot<T>(p31, u_hat);
	v3 = n_abs / length[2];
}

template<typename T>
void TriElement<T>::initialMaterials(T wf, T eps_r1, T mu_r1, T cond1, T losstangent1, T eps_r2, T mu_r2, T cond2, T losstangent2) {
	eps_out = complex<T>(eps_r1*EPS0<T>, -cond1 / wf);
	mu_out= complex<T>(mu_r1*MU0<T>, -losstangent1*mu_r1*MU0<T>);
	eps_in = complex<T>(eps_r2*EPS0<T>, -cond2 / wf);
	mu_in = complex<T>(mu_r2*MU0<T>, -losstangent2 * mu_r2*MU0<T>);
}


template<typename T>
void TriElement<T>::initialGaussPoints(int indegree, int outdegree) {
	iGPt = LookUpGaussPoints<T>(nodes, indegree);
	if (indegree != outdegree) {
		oGPt = LookUpGaussPoints<T>(nodes, outdegree);
	}
	else {
		oGPt = iGPt;
	}
}

template<typename T>
vector<T> TriElement<T>::graglia_1(int freevertex, const vector<T>& observePoint) {
	vector<T> rho = myVecMinus(observePoint, nodes[0]);
	T u0 = myVecDot(u_hat, rho);
	T v0 = myVecDot(v_hat, rho);
	T w0 = myVecDot(w_hat, rho);
	Matrix<T, 3, 1> s_minus;
	s_minus << -((length[2] - u0)*(length[2] - u3) + v0 * v3) / length[0], \
		-(u3*(u3-u0)+v3*(v3-v0))/length[1], \
		-u0;
	Matrix<T, 3, 1> s_plus;
	s_plus << ((u3 - u0)*(u3-length[2]) + v3 *(v3-v0)) / length[0], \
		 (u0*u3+v0*v3) / length[1], \
		length[2]- u0;
	Matrix<T, 3, 1> t0;
	t0 << (v0*(u3-length[2])+v3*(length[2]-u0)) / length[0], \
		(u0*v3 - v0 * u3) / length[1], \
		v0;
	Matrix<T, 3, 1> R0=(t0.array().square().array()+w0*w0).array().sqrt();
	Matrix<T, 3, 1> R_minus;
	R_minus << myVecNorm(myVecMinus(observePoint, nodes[1])), \
		myVecNorm(myVecMinus(observePoint, nodes[2])), \
		myVecNorm(myVecMinus(observePoint, nodes[0]));
	Matrix<T, 3, 1> R_plus;
	R_plus << myVecNorm(myVecMinus(observePoint, nodes[2])), \
		myVecNorm(myVecMinus(observePoint, nodes[0])), \
		myVecNorm(myVecMinus(observePoint, nodes[1]));
	Matrix<T, 3, 1> f2=((R_plus+s_plus).array()/(R_minus + s_minus).array()).log();
	if (R_plus(0) + s_plus(0) <1e-15 || R_minus(0) + s_minus(0) < 1e-15) {
		f2(0) = 0;
	}
	if (R_plus(1) + s_plus(1) < 1e-15 || R_minus(1) + s_minus(1) < 1e-15) {
		f2(1) = 0;
	}
	if (R_plus(2) + s_plus(2) < 1e-15 || R_minus(2) + s_minus(2) < 1e-15) {
		f2(2) = 0;
	}
	Matrix<T, 3, 1> f3 = (s_plus.array()*R_plus.array() - s_minus.array()*R_minus.array()) + (R0.array().square()).array()*f2.array();
	/*Matrix<T, 3, 1> beta = ((t0.array()*s_plus.array()) / (R0.array().square() + (abs(w0)*R_plus).array())).atan() - \
		((t0.array()*s_minus.array()) / (R0.array().square() + (abs(w0)*R_minus).array())).atan();*/
	//T Beta = beta.sum();
	T Beta = 0;
	if (R0(0)*R0(0) + abs(w0)*R_plus(0) != 0 && R0(0)*R0(0) + abs(w0)*R_minus(0) != 0) {
		Beta += (atan((t0(0)*s_plus(0)) / (R0(0)*R0(0) + abs(w0)*R_plus(0))) - atan((t0(0)*s_minus(0)) / (R0(0)*R0(0) + abs(w0)*R_minus(0))));
	}
	if (R0(1)*R0(1) + abs(w0)*R_plus(1) != 0 && R0(1)*R0(1) + abs(w0)*R_minus(1) != 0) {
		Beta += (atan((t0(1)*s_plus(1)) / (R0(1)*R0(1) + abs(w0)*R_plus(1))) - atan((t0(1)*s_minus(1)) / (R0(1)*R0(1) + abs(w0)*R_minus(1))));
	}
	if (R0(2)*R0(2) + abs(w0)*R_plus(2) != 0 && R0(2)*R0(2) + abs(w0)*R_minus(2) != 0) {
		Beta += (atan((t0(2)*s_plus(2)) / (R0(2)*R0(2) + abs(w0)*R_plus(2))) - atan((t0(2)*s_minus(2)) / (R0(2)*R0(2) + abs(w0)*R_minus(2))));
	}
	vector<T> ans;
	T I1 = -abs(w0)*Beta + (t0.array()*f2.array()).sum();
	T Iu = u0 * I1 + 0.5*myVecDot(u_hat, myVecPlus(myVecPlus(myVecNumProd(m_hat[0], f3[0]), myVecNumProd(m_hat[1], f3[1])), myVecNumProd(m_hat[2], f3[2])));
	T Iv = v0 * I1 + 0.5*myVecDot(v_hat, myVecPlus(myVecPlus(myVecNumProd(m_hat[0], f3[0]), myVecNumProd(m_hat[1], f3[1])), myVecNumProd(m_hat[2], f3[2])));
	if (freevertex == 0) {
		ans = myVecPlus(myVecNumProd(u_hat, Iu), myVecNumProd(v_hat, Iv));
	}
	else if (freevertex == 1) {
		ans= myVecPlus(myVecNumProd(u_hat, Iu-I1*length[2]), myVecNumProd(v_hat, Iv));
	}
	else {
		ans = myVecPlus(myVecNumProd(u_hat, Iu - I1 * u3), myVecNumProd(v_hat, Iv-I1*v3));
	}
	ans.push_back(I1);
	ans.push_back(0);
	ans.push_back(0);
	ans.push_back(0);
	if (w0 < 0) {
		Beta = -Beta;
	}
	vector<T> mf2 = myVecPlus(myVecPlus(myVecNumProd(m_hat[0],f2(0)), myVecNumProd(m_hat[1], f2(1))), myVecNumProd(m_hat[2], f2(2)));
	ans[4] = -Beta*w_hat[0]-mf2[0];
	ans[5] = -Beta * w_hat[1] - mf2[1];
	ans[6] = -Beta * w_hat[2] - mf2[2];
	return ans; //first 3 elements represent the vector of integral (rho/R), the 4th element is the integral of (1/R), the last 3 elements represent the vector of integral gradient(1/R)
	//according to Pasi Yla Oijala-2003, TAP.2003.814745, integral of rho cross gradient(1/R) = (p-q) cross integral of gradient(1/R), where p, q are free vertex of test and basis triangles.
}

template<typename T>
void RWGbasis<T>::findSharedEdges(const Mesh<T>* theMesh, const InputInfo<T> * param) {
	num_elem = theMesh->num_elem;
	triangles.resize(num_elem);
	singularThreshold = param->singularThreshold;
	wf = 2.0 * param->freq*PI<T>;
	AvgEdgeLength = 0;
	for (int i = 0; i < num_elem; ++i) {
		vector<int> M = theMesh->elems[i];
		int materialID= theMesh->attribs[i]-1; //Assuming the attribs is 1 based
		triangles[i] = new TriElement<T>(theMesh->nodes[M[0] - 1], theMesh->nodes[M[1] - 1], theMesh->nodes[M[2] - 1]);
		triangles[i]->initialMaterials(wf, param->eps_r_out[materialID], param->mu_r_out[materialID], param->cond_out[materialID], param->losstangent_out[materialID],\
			param->eps_r_in[materialID], param->mu_r_in[materialID], param->cond_in[materialID], param->losstangent_in[materialID]);
	}
	int count = 0;
	for (int i = 0; i < num_elem-1; ++i) {
		vector<int> M = theMesh->elems[i];
		for (int j = i + 1; j < num_elem; ++j) {
			vector<int> N = theMesh->elems[j];
			vector<int> tmp;
			if (M[0] == N[0] || M[1] == N[0] || M[2] == N[0]) {
				tmp.push_back(N[0]);
			}
			if (M[0] == N[1] || M[1] == N[1] || M[2] == N[1]) {
				tmp.push_back(N[1]);
			}
			if (M[0] == N[2] || M[1] == N[2] || M[2] == N[2]) {
				tmp.push_back(N[2]);
			}
			if (tmp.size() == 2) {
				if (tmp[0] != N[0]) {
					triangles[j]->edgeId[0] = count;
					triangles[j]->edgeSign[0] = -1;
				}
				else if (tmp[1] == N[1]) {
					triangles[j]->edgeId[2] = count;
					triangles[j]->edgeSign[2] = -1;
				}
				else {
					triangles[j]->edgeId[1] = count;
					triangles[j]->edgeSign[1] = -1;
				}
				if (tmp[0] != M[0] && tmp[1] != M[0]) {
					triangles[i]->edgeId[0] = count;
					triangles[i]->edgeSign[0] = 1;
				}
				else if (tmp[0] != M[1] && tmp[1] != M[1]) {
					triangles[i]->edgeId[1] = count;
					triangles[i]->edgeSign[1] = 1;
				}
				else {
					triangles[i]->edgeId[2] = count;
					triangles[i]->edgeSign[2] = 1;
				}
				count++;
				AvgEdgeLength += myVecDistance(theMesh->nodes[tmp[0]-1], theMesh->nodes[tmp[1] - 1]);
			}
		}
		
	}
	num_edge = count;
	AvgEdgeLength /= (T)num_edge;
}

template<typename T>
void RWGbasis<T>::GaussPoints(int indegree, int outdegree) {
	innerGaussWeight = LookUpGaussWeights<T>(indegree);
	outerGaussWeight = LookUpGaussWeights<T>(outdegree);
	for (int i = 0; i < num_elem; ++i) {
		TriElement<T>* tmp = triangles[i];
		tmp->initialGaussPoints(indegree, outdegree);
	}
}



#endif
