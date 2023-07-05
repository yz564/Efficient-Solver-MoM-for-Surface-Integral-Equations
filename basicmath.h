#ifndef BASICMATH_H// include guard
#define BASICMATH_H
#include<assert.h>
#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

template<typename T>
const T PI = 3.141592653589793238463;
template<typename T>
const T EPS0 = 8.8541878128131313e-12;
template<typename T>
const T MU0 = 1.25663706212191919e-6;

template<typename T>
vector<T> myVecPlus(const vector<T> & x, const vector<T> & y) {
	assert(x.size() == y.size());
	vector<T> ans(x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		ans[i]=(x[i]+y[i]);
	}
	return ans;
}

template<typename T>
vector<T> myVecMinus(const vector<T> & x, const vector<T> & y) {
	assert(x.size() == y.size());
	vector<T> ans(x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		ans[i]=(x[i] - y[i]);
	}
	return ans;
}

template<typename T>
vector<complex<T> > myVecReal2Complex(const vector<T> & x) {
	vector<complex<T> > ans;
	for (size_t i = 0; i < x.size(); ++i) {
		ans.push_back(complex<T>(x[i], 0));
	}
	return ans;
}

template<typename T>
Matrix <T, Dynamic, Dynamic> myVec2Eigen(const vector<T > & x) {
	Matrix <T, Dynamic, Dynamic> ans(1, x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		ans(0,i) = x[i];
	}
	return ans;
}

template<typename T>
complex<T> myReal2Complex(const T & x) {
	complex<T> ans= complex<T>(x,0);
	return ans;
}

template<typename T>
vector<T> myVecNumProd(const vector<T> & x, const T a) {
	vector<T> ans(x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		ans[i]=x[i]*a;
	}
	return ans;
}

template<typename T>
vector<T> myVecNumDiv(const vector<T> & x, const T a) {
	vector<T> ans(x.size());
	for (size_t i = 0; i < x.size(); ++i) {
		ans[i] = x[i]/a;
	}
	return ans;
}


template<typename T>
T myVecDot(const vector<T> & x, const vector<T> & y) {
	assert(x.size() == y.size());
	T ans = 0;
	for (size_t i = 0; i < x.size(); ++i) {
		ans += x[i] * y[i];
	}
	return ans;
}

template<typename T>
vector<T> myVecCross(const vector<T> & x, const vector<T> & y) {
	assert(x.size() == y.size() && x.size()==3);
	vector<T> ans(3);
	ans[0] = x[1] * y[2] - x[2] * y[1];
	ans[1] = x[2] * y[0] - x[0] * y[2];
	ans[2] = x[0] * y[1] - x[1] * y[0];
	return ans;
}


template<typename T>
T myVecNorm(const vector<T> & x) {
	T ans=0;
	for (size_t i = 0; i < x.size(); ++i) {
		ans += x[i] * x[i];
	}
	return sqrt(ans);
}

template<typename T>
T myVecDistance(const vector<T> & x, const vector<T> & y) {
	T ans = 0;
	T tmp = 0;
	for (size_t i = 0; i < x.size(); ++i) {
		tmp = x[i] - y[i];
		ans += (tmp*tmp);
	}
	return sqrt(ans);
}

template<typename T>
T myVecMean(const vector<T> & x) {
	T ans=0;
	if (x.size() == 0) {
		return ans;
	}
	for (size_t i = 0; i < x.size(); ++i) {
		ans += x[i];
	}
	return ans/x.size();
}

template<typename T>
T myVecMax(const vector<T> & x) {
	if (x.size() == 0) {
		return 0;
	}
	T ans = x[0];
	for (size_t i = 1; i < x.size(); ++i) {
		if (ans < x[i]) {
			ans = x[i];
		}
	}
	return ans;
}

template<typename T>
T myVecMin(const vector<T> & x) {
	if (x.size() == 0) {
		return 0;
	}
	T ans = x[0];
	for (size_t i = 1; i < x.size(); ++i) {
		if (ans > x[i]) {
			ans = x[i];
		}
	}
	return ans;
}


#endif