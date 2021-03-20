//Program for adding some functions for Eigen
//(Functions like mat_ind is unnecessary in C++ because loop calculation is efficient in C++)
///////////////////////////////////////////////////////
//By Wenaho Xu, Xi'an Jiaotong University
#ifndef EIGEN_PLUS_V1_H
#define EIGEN_PLUS_V1_H
#define EIGEN_USE_MKL_ALL

#include <algorithm>
#include <complex>
#include "cpp_plus_v1.h"
#include <cstdlib>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <unsupported/Eigen/IterativeSolvers>

using namespace std;
using namespace Eigen;

//Class for Eigen matrix-free matrix-vector multiplication
//class MatrixReplacement;
//namespace Eigen {
//	namespace internal{
//		template<>
//		//(The following declaration doesn't affect global matrix and right-hand-side terms)
//		struct traits<MatrixReplacement> : public internal::traits<SparseMatrix<double> > {};
//	}
//}
//class MatrixReplacement : public EigenBase<MatrixReplacement>
//{
//public:
//	typedef double Scalar;
//	typedef double RealScalar;
//	typedef int StorageIndex;
//	enum
//	{
//		ColsAtCompileTime = Dynamic,
//		MaxColsAtCompileTime = Dynamic,
//		IsRowMajor = false //This parameter doesn't affect global matrix
//	};
//	int row_num, col_num;
//public:
//	MatrixReplacement(int row_num, int col_num)
//	{
//		this->row_num = row_num;
//		this->col_num = col_num;
//	}
//	Index rows() const { return row_num; }
//	Index cols() const { return col_num; }
//	template<typename Rhs>
//	Product<MatrixReplacement, Rhs, AliasFreeProduct> operator*
//		(const MatrixBase<Rhs>& x) const
//	{
//		return Product<MatrixReplacement, Rhs, AliasFreeProduct>(*this, x.derived());
//	}
//};

//Add the function of general clockwise process
template <typename T>
Matrix<T, -1, -1> cwise_process(Matrix<T, -1, -1> A, T process_fun(T))
{
	Matrix<T, -1, -1> B = A;
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < A.size(); i++)
	{
		B(i) = process_fun(A(i));
	}
	return B;
}

//Add function for clockwise pow_int
template <typename T>
T cwise_pow_int(T A, int n)
{
	T B = A;
	int N = (int)A.size();
	int i;
#pragma omp parallel for private (i)
	for (i = 0; i < N; i++)
	{
		B(i) = pow_int(A(i), n);
	}
	return B;
}

//Function for getting coordinate (COO) format from Eigen sparse matrix
template <typename T>
void sp2coo(SparseMatrix<T, RowMajor> &A, RowVectorXi &row_ind, RowVectorXi &col_ind,
	Matrix<T, 1, -1> &val)
{
	int row_num = (int)A.rows();
	int nnz = (int)A.nonZeros();
	RowVectorXi row_first_ind(row_num + 1);
	row_ind = RowVectorXi(nnz);
	col_ind = RowVectorXi(nnz);
	val = Matrix<T, 1, -1>(nnz);
	memcpy(&row_first_ind(0), A.outerIndexPtr(), (row_num + 1) * sizeof(int));
	memcpy(&col_ind(0), A.innerIndexPtr(), nnz * sizeof(int));
	memcpy(&val(0), A.valuePtr(), nnz * sizeof(T));
	int i, row_nonzero_num;
	int count = 0;
	for (i = 0; i < row_num; i++)
	{
		row_nonzero_num = row_first_ind(i + 1) - row_first_ind(i);
		if (row_nonzero_num != 0)
		{
			row_ind.segment(count, row_nonzero_num).setConstant(i);
			count += row_nonzero_num;
		}
	}
}

//Function for getting Eigen sparse matrix by coordinate (COO) format
template<typename T>
SparseMatrix<T, RowMajor> coo2sp(RowVectorXi row_ind, RowVectorXi col_ind, Matrix<T, 1, -1> &val,
	int row_num, int col_num)
{
	int k;
	int nnz = (int)val.size();
	vector<Triplet<T> > triplet_list(nnz);
#pragma omp parallel for private(k)
	for (k = 0; k < nnz; k++)
	{
		triplet_list[k] = Triplet<T>(row_ind(k), col_ind(k), val(k));
	}
	SparseMatrix<T, RowMajor> A = SparseMatrix<T, RowMajor>(row_num, col_num);
	A.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return A;
}

//Add the funciton of 'amd' (approximate minimum degree)
template <typename T>
VectorXi amd(SparseMatrix<T, RowMajor> &A)
{
	SparseMatrix<T, ColMajor> ColMajor_A = A;
	AMDOrdering<int> ordering;
	PermutationMatrix<Dynamic, Dynamic, int> perm;
	ordering(ColMajor_A, perm);
	return perm.indices();
}

//Function for getting the inverse ordering vector
VectorXi get_inv_ordering_vec(VectorXi ordering_vec)
{
	int n = (int)ordering_vec.size();
	VectorXi inv_ordering_vec(n);
	int i;
	for (i = 0; i < n; i++)
	{
		inv_ordering_vec((ordering_vec(i))) = i;
	}
	return inv_ordering_vec;
}

//Function for clockwise conversion
template <typename T1, typename T2>
Matrix<T1, -1, -1> cwise_convert(T2 A)
{
	int row_num = (int)A.rows();
	int col_num = (int)A.cols();
	Matrix<T1, -1, -1> A1(row_num, col_num);
	int i, j;
#pragma omp parallel for private(j,i)
	for (j = 0; j < row_num; j++)
	{
		for (i = 0; i < col_num; i++)
		{
			A1(j, i) = (T1)A(j, i);
		}
	}
	return A1;
}
template <typename T1, typename T2>
SparseMatrix<T1, RowMajor> cwise_convert(SparseMatrix<T2, RowMajor> A)
{
	int row_num = (int)A.rows();
	int col_num = (int)A.cols();
	RowVectorXi row_ind, col_ind;
	Matrix<T2, 1, -1> val;
	sp2coo<T2>(A, row_ind, col_ind, val);
	int nnz = (int)A.nonZeros();
	Matrix<T1, 1, -1> converted_val(nnz);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < nnz; i++)
	{
		converted_val(i) = (T1)val(i);
	}
	return coo2sp<T1>(row_ind, col_ind, converted_val, row_num, col_num);
}

//Function for geometric row scaling
template <typename T>
void geometric_row_scaling(T *A_ptr, VectorXd &A_row_norm_vec)
{
	int row_num = (int)A_ptr->rows();
	A_row_norm_vec = VectorXd(row_num);
	int j;
#pragma omp parallel for private(j)
	for (j = 0; j < row_num; j++)
	{
		A_row_norm_vec(j) = A_ptr->row(j).norm();
		A_ptr->row(j) /= A_row_norm_vec(j);
	}
}

//Function for getting maximum amplitdue of sparse matrix
template <typename T>
double get_max_amp(SparseMatrix<T, RowMajor> A)
{
	int nnz = A.nonZeros();
	Matrix<T, 1, -1> val(nnz);
	memcpy(&val(0), A.valuePtr(), nnz * sizeof(T));
	return val.cwiseAbs().maxCoeff();
}

//Function for getting average amplitdue of sparse matrix
template <typename T>
double get_ave_nonzero_amp(SparseMatrix<T, RowMajor> A)
{
	int nnz = A.nonZeros();
	Matrix<T, 1, -1> val(nnz);
	memcpy(&val(0), A.valuePtr(), nnz * sizeof(T));
	return val.cwiseAbs().mean();
}

//Function for reducing nonzero elements in sparse matrix
template<typename T>
SparseMatrix<T, RowMajor> reduce_nonzero_elem(SparseMatrix<T, RowMajor> A, double elem_reduce_ratio)
{
	int nnz = (int)A.nonZeros();
	VectorXd row_scaling = get_row_norm(A).cwiseInverse();
	RowVectorXi row_ind, col_ind;
	Matrix<T, 1, -1> val, backup_val;
	sp2coo(A, row_ind, col_ind, backup_val);
	apply_row_scaling(A, row_scaling);
	sp2coo(A, row_ind, col_ind, val);
	RowVectorXd abs_val(nnz), backup_abs_val;
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < nnz; i++)
	{
		abs_val(i) = (double)abs(val(i));
	}
	backup_abs_val = abs_val;
	sort_by_key(&abs_val(0), &row_ind(0), nnz, false);
	abs_val = backup_abs_val;
	sort_by_key(&abs_val(0), &col_ind(0), nnz, false);
	abs_val = backup_abs_val;
	val = backup_val;
	sort_by_key(&abs_val(0), &val(0), nnz, false);
	int left_elem_num = (int)floor(nnz*(1 - elem_reduce_ratio));
	return coo2sp<T>((RowVectorXi)row_ind.head(left_elem_num), (RowVectorXi)col_ind.head(left_elem_num),
		(Matrix<T, 1, -1>)val.head(left_elem_num), (int)A.rows(), (int)A.cols());
}

//Add the function of 'unique'
template <typename T>
Matrix<T, -1, 1> unique(Matrix<T, -1, -1> A)
{
	int L = (int)A.size();
	sort(&A(0), &A(0) + L);
	T *p = unique(&A(0), &A(0) + L);
	int len = p - &A(0);
	Matrix<T, -1, 1> a(len);
	memcpy(&a(0), &A(0), len * sizeof(T));
	return a;
}

//Add the function of 'cat(dim,A,B)'
template <typename T>
T cat(int dim, T A, T B)
{
	int m1 = (int)(A.rows()), n1 = (int)(A.cols());
	int m2 = (int)(B.rows()), n2 = (int)(B.cols());
	T C;
	if (dim == 1) //Concatinate in the row direction
	{
		int i;
		if (n1 != n2)
		{
			printf("\nThe number of column must be same\n");
			exit(1);
		}
		C = T(m1 + m2, n1);
		for (i = 0; i < m1; i++)
		{
			C.row(i) = A.row(i);
		}
		for (i = 0; i < m2; i++)
		{
			C.row(m1 + i) = B.row(i);
		}
	}
	else if (dim == 2) //Concatinate in the column direction
	{
		int j;
		if (m1 != m2)
		{
			printf("\nThe number of row must be same\n");
			exit(1);
		}
		C = T(m1, n1 + n2);
		for (j = 0; j < n1; j++)
		{
			C.col(j) = A.col(j);
		}
		for (j = 0; j < n2; j++)
		{
			C.col(n1 + j) = B.col(j);
		}
	}
	else
	{
		printf("\nThe direction of the concatenation is incorrect\n\n");
		exit(1);
	}
	return C;
}

//Add the function of 'repmat'
template <typename T>
T repmat(T A, int m, int n)
{
	//Initialize extended matrix
	int m0 = (int)(A.rows());
	int n0 = (int)(A.cols());
	T big_A = T::Zero(m0*m, n0*n);

	//Pad extended matrix with blocks
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			big_A.block(i*m0, j*n0, m0, n0) = A;
		}
	}

	return big_A;
}

//Add the function of 'mod'
MatrixXi mod(MatrixXi A, int N)
{
	int i;
	for (i = 0; i < (int)A.size(); i++)
	{
		A(i) %= N;
		A(i) = (A(i) + N) % N;
	}
	return A;
}
MatrixXd mod(MatrixXd A, double x)
{
	int i;
	for (i = 0; i < (int)A.size(); i++)
	{
		A(i) = fmod(A(i), x);
		A(i) = fmod(A(i) + x, x);
	}
	return A;
}

//Add the function of 'fftshfit(A,dim)'
template <typename T>
T fftshift(T A, int dim)
{
	int m = (int)(A.rows());
	int n = (int)(A.cols());
	T B(m, n);
	if (dim == 1)
	{
		B.topRows(m / 2) = A.bottomRows(m / 2);
		B.bottomRows(m - m / 2) = A.topRows(m - m / 2);
	}
	else if (dim == 2)
	{
		B.leftCols(n / 2) = A.rightCols(n / 2);
		B.rightCols(n - n / 2) = A.leftCols(n - n / 2);
	}
	else
	{
		printf("\nThe dim of fftshift is incorrect\n");
		exit(1);
	}
	return B;
}

//Add the function of 'ifftshfit(A,dim)'
template <typename T>
T ifftshift(T B, int dim)
{
	int m = (int)(B.rows());
	int n = (int)(B.cols());
	T A(m, n);
	if (dim == 1)
	{
		A.bottomRows(m / 2) = B.topRows(m / 2);
		A.topRows(m - m / 2) = B.bottomRows(m - m / 2);
	}
	else if (dim == 2)
	{
		A.rightCols(n / 2) = B.leftCols(n / 2);
		A.leftCols(n - n / 2) = B.rightCols(n - n / 2);
	}
	else
	{
		printf("\nThe dim of ifftshift is incorrect\n");
		exit(1);
	}
	return A;
}

//Add the function of 'circshift(A,K,dim)'
template<typename T>
T circshift(T A, int K, int dim = 1)
{
	int row_num = (int)(A.rows());
	int col_num = (int)(A.cols());
	T B(row_num, col_num);
	if (dim == 1)
	{
		if (K > 0)
		{
			B.bottomRows(row_num - K) = A.topRows(row_num - K);
			B.topRows(K) = A.bottomRows(K);
		}
		else
		{
			B.topRows(row_num + K) = A.bottomRows(row_num + K);
			B.bottomRows(-K) = A.topRows(-K);
		}
	}
	else if (dim == 2)
	{
		if (K > 0)
		{
			B.rightCols(col_num - K) = A.leftCols(col_num - K);
			B.leftCols(K) = A.rightCols(K);
		}
		else
		{
			B.leftCols(col_num + K) = A.rightCols(col_num + K);
			B.rightCols(-K) = A.leftCols(-K);
		}
	}
	else
	{
		printf("\nThe dim of circshift is incorrect\n");
		exit(1);
	}
	return B;
}

//Add the funciton of getting MP inverse, using SVD decomposition
template <typename T>
T svd_MP_inv(T A)
{
	int m = (int)(A.rows());
	JacobiSVD<T> svd(A, ComputeThinU | ComputeThinV);
	return svd.solve(T::Identity(m, m));
}

//Add the function of 'sub2ind', where the objective matrix is column-major
// (Notice that row_sub, col_sub and returned indices all start from 0 
//  so as to be consistent with the setting in C++)
template <typename T>
T sub2ind(int row_num, int col_num, T row_sub, T col_sub)
{

	//Check input
	if (row_sub.size() != col_sub.size())
	{
		cout << "The row_sub and col_sub must be in the same size!" << endl;
		exit(1);
	}
	if (row_sub.maxCoeff() >= row_num)
	{
		cout << "The element of row_sub must be smaller than row_num!" << endl;
		exit(1);
	}
	if (col_sub.maxCoeff() >= col_num)
	{
		cout << "The element of col_sub must be smaller than col_num" << endl;
		exit(1);
	}

	//Get required indices
	int m = (int)(row_sub.rows());
	int n = (int)(row_sub.cols());
	T ind_mat(m, n);
	int i; int j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			ind_mat(i, j) = col_sub(i, j)*row_num + row_sub(i, j);
		}
	}

	return ind_mat;
}

//Add the function of 'fliplr'
template <typename T>
T fliplr(T A)
{
	return A.rowwise().reverse();
}

//Add the function of 'flipud'
template <typename T>
T flipud(T A)
{
	return A.colwise().reverse();
}

// Add the function of chasing method
///////////////////////////////////////////////////////
//a is the -1-th diagonal of the tridiagonal coefficient matrix
//b is the 0-th (main) diagonal of the tridiagonal coefficient matrix
//c is the 1-th diagonal of the tridiagonal coefficient matrix
//d is the right-hand vector
///////////////////////////////////////////////////////
void tridiagonal_LU(VectorXd a, VectorXd b, VectorXd c, VectorXd &iota, VectorXd &u)
{
	//Check input
	int n = (int)b.size();
	if (a.size() != n - 1 || c.size() != n - 1)
	{
		cout << "tridiagonal_LU: " << endl;
		cout << "the sizes of of a and c must be equal to n-1!" << endl;
		exit(1);
	}

	//Execute LU decomposition
	iota = VectorXd(n - 1);
	u = VectorXd(n);
	u(0) = b(0);
	int i, k;
	for (i = 1; i < n; i++)
	{
		k = i - 1;
		iota(k) = a(k) / u(i - 1);
		u(i) = b(i) - iota(k)*c(i - 1);
	}
}
template <typename T>
T tridiagonal_chasing(VectorXd &iota, VectorXd &u, VectorXd &c, T d)
{
	//Check input
	int n = (int)u.size();
	if (iota.size() != n - 1 || c.size() != n - 1 || d.size() != n)
	{
		cout << "tridiagonal_chasing: " << endl;
		cout << "The sizes of u and d must be equal to n, and the sizes of of iota and c must \
be equal to n-1!" << endl;
		exit(1);
	}

	//Solving linear equations with chasing method
	T y(n), x(n);
	y(0) = d(0);
	int i, k;
	for (i = 1; i < n; i++)
	{
		k = i - 1;
		y(i) = d(i) - iota(k)*y(i - 1);
	}
	x(n - 1) = y(n - 1) / u(n - 1);
	for (i = n - 2; i > -1; i--)
	{
		x(i) = (y(i) - c(i)*x(i + 1)) / u(i);
	}

	return x;
}
template<typename T>
T chasing_method(VectorXd a, VectorXd b, VectorXd c, T d)
{
	T iota, u;
	tridiagonal_LU(a, b, c, iota, u);
	return tridiagonal_chasing(iota, u, c, d);
}

//Add the function of 'interp1' with 'linear' option, where the original x is assumed
// to start from x0 with uniform interval h
template <typename T>
T simple_linear_interp1(double x0, double h, T y, VectorXd xi)
{
	//Check input
	int n = (int)y.size();
	int n1 = (int)xi.size();
	if ((int)y.rows() != 1 && (int)y.cols() != 1)
	{
		cout << "y must be a vector rather than a matrix!" << endl;
		exit(1);
	}

	//Initialization
	T yi(n1);
	int k;
	RowVectorXi i_vec(n1);

	//Execute 1D linear interpolation
#pragma omp parallel for private(k)
	for (k = 0; k < n1; k++)
	{
		i_vec(k) = min(max((int)ceil((xi(k) - x0) / h), 1), n - 1);
		yi(k) = (x0 + i_vec(k)*h - xi(k)) / h * y(i_vec(k) - 1) +
			(xi(k) - x0 - (i_vec(k) - 1)*h) / h * y(i_vec(k));
	}

	return yi;
}

//Add the function of 'interp2' with 'linear' option, where the original x and y are assumed
// to start from x0 and y0 with uniform interval dx and dy, respectively
template <typename T>
T simple_linear_interp2(double x0, double y0, double dx, double dy, Matrix<T, -1, -1> &v,
	double xi, double yi)
{
	//Check input
	int nx = (int)v.cols();
	int ny = (int)v.rows();
	if (nx == 1 || ny == 1)
	{
		cout << "v must be a matrix rather than a vector!" << endl;
		exit(1);
	}

	//Initialization
	int i, j;

	//Execute 2D linear interpolation
	i = (int)ceil((xi - x0) / dx);
	i = max(i, 1);
	i = min(i, nx - 1);
	j = (int)ceil((yi - y0) / dy);
	j = max(j, 1);
	j = min(j, ny - 1);
	// (Do 1D interpolation in the y-direction first, and then do 1D interpolation in the x-direction)
	return (x0 + i * dx - xi) / dx * ((y0 + j * dy - yi) / dy * v(j - 1, i - 1) +
		(yi - y0 - (j - 1)*dy) / dy * v(j, i - 1))
		+ (xi - x0 - (i - 1)*dx) / dx * ((y0 + j * dy - yi) / dy * v(j - 1, i) +
		(yi - y0 - (j - 1)*dy) / dy * v(j, i));
}

//Add function for 3D resampling by linear interpolation
template <typename T>
Matrix<Matrix<T, -1, -1>, -1, 1> linear_3D_resample(double dx, double dy, double dz,
	Matrix<Matrix<T, -1, -1>, -1, 1> v,
	double x1_start, double y1_start, double z1_start,
	double dx1, double dy1, double dz1,
	int nx1, int ny1, int nz1)
{
	//Initialization
	int nx = (int)v(0).cols();
	int ny = (int)v(0).rows();
	int nz = (int)v.size();
	Matrix<Matrix<T, -1, -1>, -1, 1> tmp_v;
	int i, j, k;

	//Do interpolation in the z direction
	VectorXd z1 = z1_start + VectorXd::LinSpaced(nz1, 0, nz1 - 1).array()*dz1;
	tmp_v = simple_linear_interp1(0, dz, v, z1);

	//Do interpolation in the x- and y- directions
	Matrix<Matrix<T, -1, -1>, -1, 1> v1(nz1);
	for (k = 0; k < nz1; k++)
	{
		v1(k) = Matrix<T, -1, -1>(ny1, nx1);
#pragma omp parallel for private(j,i)
		for (j = 0; j < ny1; j++)
		{
			for (i = 0; i < nx1; i++)
			{
				v1(k)(j, i) = simple_linear_interp2(0, 0, dx, dy, tmp_v(k),
					x1_start + i * dx1, y1_start + j * dy1);
			}
		}
	}

	return v1;
}
template <typename T>
Matrix<Matrix<T, -1, -1>, -1, 1> linear_3D_resample(double dx, double dy, double dz,
	Matrix<Matrix<T, -1, -1>, -1, 1> v,
	double x1_start, double y1_start, double z1_start,
	double dx1, double dy1, double dz1, double int_fun(double))
{
	//Initialization
	int nx = (int)v(0).cols();
	int ny = (int)v(0).rows();
	int nz = (int)v.size();
	int nx1 = (int)int_fun((nx - 1)*dx / dx1 + 1);
	int ny1 = (int)int_fun((ny - 1)*dy / dy1 + 1);
	int nz1 = (int)int_fun((nz - 1)*dz / dz1 + 1);
	return linear_3D_resample(dx, dy, dz, v, x1_start, y1_start, z1_start, dx1, dy1, dz1, nx1, ny1, nz1);
}

//Add the function of 'interp1' with 'nearest' option, where the original x is assumed
// to start from x0 with uniform interval h
template <typename T>
T simple_nearest_interp1(double x0, double h, T y, VectorXd xi)
{
	//Check input
	int n = (int)y.size();
	int n1 = (int)xi.size();
	if ((int)y.rows() != 1 && (int)y.cols() != 1)
	{
		cout << "y must be a vector rather than a matrix!" << endl;
		exit(1);
	}

	//Initialization
	T yi(n1);
	double x;
	int i, k;

	//Execute 1D linear interpolation
	for (k = 0; k < n1; k++)
	{
		i = (int)round((xi(k) - x0) / h);
		i = max(i, 0);
		i = min(i, n - 1);
		yi(k) = y(i);
	}

	return yi;
}

//Add the function of 'randn(row_num,col_num)'
MatrixXd randn(int row_num, int col_num)
{
	//Initialization
	srand((unsigned int)time(0));
	MatrixXd A(row_num, col_num);
	MatrixXd U1 = MatrixXd::Random(row_num, col_num);
	MatrixXd U2 = MatrixXd::Random(row_num, col_num);
	double epsilon = DBL_EPSILON;
	double counter_epsilon = 1 - epsilon;
	double pi = acos(-1.0);
	int i;

	//Get normal distribution by "Box-Muller" method
	U1 = (U1.array() + 1) / 2;
	for (i = 0; i < U1.size(); i++)
	{
		if (U1(i) < epsilon)
		{
			U1(i) = epsilon;
		}
		if (U1(i) > counter_epsilon)
		{
			U1(i) = counter_epsilon;
		}
	}
	U2 = (U2.array() + 1) / 2;
	A = (-2 * U1.array().log()).sqrt()*(2 * pi*U2).array().cos();

	return A;
}

//Add the function of awgn(A,SNR,'measured')
MatrixXd awgn(MatrixXd A, double SNR)
{

	//Get the standard deviation of Gaussian white noise
	double signal_power = A.squaredNorm() / A.size();
	double noise_power = signal_power / pow(10.0, SNR / 10);
	double sigma = sqrt(noise_power);

	//Add noise
	MatrixXd noise = sigma * randn((int)A.rows(), (int)A.cols());

	return A + noise;
}

//Add the function for analytically solving Vandermonde linear equations V(x0,...,xn)*z=b
//////////////////////////////////////////////////////////////////////
//x is the vector that constructs Vandermonde coefficient matrix
//b is the right-hand side term
//////////////////////////////////////////////////////////////////////
//z is the returned solution
VectorXd Golub_solve_Van(VectorXd x, VectorXd b)
{
	//Check input
	if (x.size() != b.size())
	{
		cout << "x and b must be in the same size!" << endl;
		exit(1);
	}

	//Initialization
	int n = (int)x.size() - 1;
	VectorXd z = b;
	VectorXd z1(n + 1);
	int k;

	//Get analytic solution by recursive algorithm
	for (k = 0; k <= n - 1; k++)
	{
		memcpy(&(z1(k)), &(z(k)), (n - k) * sizeof(double));
		z.tail(n - k) -= x(k)*z1.segment(k, n - k);
	}
	for (k = n - 1; k >= 0; k--)
	{
		z.tail(n - k) = z.tail(n - k).cwiseQuotient(x.tail(n - k) - x.head(n - k));
		z.segment(k, n - k) -= z.tail(n - k);
	}

	return z;
}

//Add max, min and mean functions for high-dimensional matrix
double max(Matrix<MatrixXd, -1, 1> &A)
{
	double max_val = A(0).maxCoeff();
	int i;
	for (i = 1; i < A.size(); i++)
	{
		max_val = max(max_val, A(i).maxCoeff());
	}
	return max_val;
}
double min(Matrix<MatrixXd, -1, 1> &A)
{
	double min_val = A(0).minCoeff();
	int i;
	for (i = 1; i < A.size(); i++)
	{
		min_val = min(min_val, A(i).minCoeff());
	}
	return min_val;
}
double mean(Matrix<MatrixXd, -1, 1> &A)
{
	double sum = A(0).sum();
	int n = (int)A(0).size();
	int i;
	for (i = 1; i < A.size(); i++)
	{
		sum += A(i).sum();
		n += (int)A(i).size();
	}
	return sum / n;
}

//Add function for getting the vector in the third dimension, which I call pillar
VectorXd get_pillar(Matrix<MatrixXd, -1, 1> &A, int row_ind, int col_ind)
{
	int nz = (int)A.size();
	VectorXd pillar(nz);
	int k;
	for (k = 0; k < nz; k++)
	{
		pillar(k) = A(k)(row_ind, col_ind);
	}
	return pillar;
}

//Add function for setting the vector in the third dimension, which I call pillar
void set_pillar(Matrix<MatrixXd, -1, 1> &A, int row_ind, int col_ind, VectorXd pillar)
{
	//Check input
	int nz = (int)A.size();
	if (pillar.size() != nz)
	{
		cout << "set_pillar: " << endl;
		cout << "The size of pillar must be consistent with the size of A!" << endl;
	}

	//Set pillar
	int k;
	for (k = 0; k < nz; k++)
	{
		A(k)(row_ind, col_ind) = pillar(k);
	}
}

//Add function for transforming absolute wave impedance to reflectivity
MatrixXd imp2R(MatrixXd imp)
{
	//Initialization
	int n = (int)imp.rows();
	int trace_num = (int)imp.cols();
	MatrixXd R(n, trace_num);

	//Calculate reflectivity
	R.topRows(n - 1) = (imp.bottomRows(n - 1) - imp.topRows(n - 1)).array() /
		((imp.bottomRows(n - 1) + imp.topRows(n - 1)).array() + DBL_EPSILON);
	R.row(n - 1) = R.row(n - 2);

	return R;
}

//Add function for getting the difference of indices
template <typename T>
T ind_diff(int ind_num, T cur_ind)
{
	int i, j;
	int cur_ind_num = (int)cur_ind.size();
	int diff_ind_num = ind_num - cur_ind_num;
	T diff_ind(diff_ind_num);
	int count = 0;
	for (i = 0; i < cur_ind(0); i++)
	{
		diff_ind(count) = i;
		count++;
	}
	for (i = 0; i < cur_ind_num - 1; i++)
	{
		for (j = cur_ind(i) + 1; j < cur_ind(i + 1); j++)
		{
			diff_ind(count) = j;
			count++;
		}
	}
	for (i = cur_ind(cur_ind_num - 1) + 1; i < ind_num; i++)
	{
		diff_ind(count) = i;
		count++;
	}
	return diff_ind;
}

//Add functions for setting zero of a row or a column of a 3D matrix
void set_row_zero(Matrix<VectorXd, -1, -1> &A, int row_ind)
{
	int i;
	for (i = 0; i < A.cols(); i++)
	{
		A(row_ind, i).setZero();
	}
}
void set_col_zero(Matrix<VectorXd, -1, -1> &A, int col_ind)
{
	int j;
	for (j = 0; j < A.rows(); j++)
	{
		A(j, col_ind).setZero();
	}
}

//Add functions for copying a row or a column of a 3D matrix
void copy_row(Matrix<VectorXd, -1, -1> &A, Matrix<VectorXd, -1, -1> &B, int row_ind)
{
	int i;
	for (i = 0; i < A.cols(); i++)
	{
		A(row_ind, i) = B(row_ind, i);
	}
}
void copy_col(Matrix<VectorXd, -1, -1> &A, Matrix<VectorXd, -1, -1> &B, int col_ind)
{
	int j;
	for (j = 0; j < A.rows(); j++)
	{
		A(j, col_ind) = B(j, col_ind);
	}
}

//Add function for transformations between pillar 2D matrix and pillar 3D matrix
template <typename T>
Matrix<Matrix<T, -1, 1>, -1, -1> pillar_2D_2_pillar_3D(Matrix<T, -1, -1> S2,
	int inline_num, int xline_num)
{
	int trace_num = (int)S2.cols();
	if (trace_num != inline_num * xline_num)
	{
		cout << "pillar_2D_2_pillar_3D: trace_num must be equal to inline_num*xline_num!" << endl;
		exit(1);
	}
	Matrix<Matrix<T, -1, 1>, -1, -1> S3(xline_num, inline_num);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < trace_num; i++)
	{
		S3(i) = S2.col(i);
	}
	return S3;
}
template <typename T>
Matrix<T, -1, -1> pillar_3D_2_pillar_2D(Matrix<Matrix<T, -1, 1>, -1, -1> S3)
{
	int trace_num = (int)S3.size();
	Matrix<T, -1, -1> S2(S3(0).size(), trace_num);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < trace_num; i++)
	{
		S2.col(i) = S3(i);
	}
	return S2;
}

//Add functions for transformations between pillar 2D matrix and layered 3D matrix
template <typename T>
Matrix<Matrix<T, -1, -1>, -1, 1> pillar_2D_2_layered_3D(Matrix<T, -1, -1> S2,
	int inline_num, int xline_num)
{
	int trace_num = (int)S2.cols();
	int nt = (int)S2.rows();
	if (trace_num != inline_num * xline_num)
	{
		cout << "pillar_2D_2_layered_3D: trace_num must be equal to inline_num*xline_num!" << endl;
		exit(1);
	}
	S2.transposeInPlace();
	Matrix<Matrix<T, -1, -1>, -1, 1> S3(nt);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < nt; i++)
	{
		S3(i) = Matrix<T, -1, -1>(xline_num, inline_num);
		memcpy(&(S3(i)(0)), &S2(0, i), trace_num * sizeof(T));
	}
	S2.transposeInPlace();
	return S3;
}
template <typename T>
Matrix<T, -1, -1> layered_3D_2_pillar_2D(Matrix<Matrix<T, -1, -1>, -1, 1> S3)
{
	int nt = (int)S3.size();
	int trace_num = (int)S3(0).size();
	Matrix<T, -1, -1> S2(trace_num, nt);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < nt; i++)
	{
		memcpy(&S2(0, i), &(S3(i)(0)), trace_num * sizeof(T));
	}
	return S2.transpose();
}

//Add functions for transformations between pillar 2D matrix and line 3D matrix
template <typename T>
Matrix<Matrix<T, -1, -1>, 1, -1> pillar_2D_2_line_3D(Matrix<T, -1, -1> S2,
	int inline_num, int xline_num)
{
	int trace_num = (int)S2.cols();
	if (trace_num != inline_num * xline_num)
	{
		cout << "pillar_2D_2_line_3D: trace_num must be equal to inline_num*xline_num!"
			<< endl;
		exit(1);
	}
	Matrix<Matrix<T, -1, -1>, 1, -1> S3(inline_num);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < inline_num; i++)
	{
		S3(i) = S2.middleCols(i*xline_num, xline_num);
	}
	return S3;
}
template <typename T>
Matrix<T, -1, -1> line_3D_2_pillar_2D(Matrix<Matrix<T, -1, -1>, 1, -1> S3)
{
	int nt = (int)S3(0).rows();
	int inline_num = (int)S3.size();
	int xline_num = (int)S3(0).cols();
	Matrix<T, -1, -1> S2(nt, inline_num*xline_num);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < inline_num; i++)
	{
		S2.middleCols(i*xline_num, xline_num) = S3(i);
	}
	return S2;
}

//Add functions for transformations between pillar 1D vector and layered 3D matrix
template <typename T>
Matrix<Matrix<T, -1, -1>, -1, 1> pillar_1D_2_layered_3D(Matrix<T, -1, 1> s,
	int inline_num, int xline_num, int nt)
{
	//Check input
	if (s.size() != (long int)inline_num*(long int)xline_num*(long int)nt)
	{
		cout << "pillar_1D_2_layered_3D: the size of s must be equal to inline_num*xline_num*nt!" << endl;
		exit(1);
	}

	//Do transformation
	Matrix<Matrix<T, -1, -1>, -1, 1> S3(nt);
	int i;
	int layer_size = inline_num * xline_num;
#pragma omp parallel for private(i)
	for (i = 0; i < nt; i++)
	{
		S3(i) = Matrix<T, -1, -1>(xline_num, inline_num);
		memcpy(&(S3(i)(0)), &s(0 + i * layer_size), layer_size * sizeof(T));
	}

	return S3;
}
template <typename T>
Matrix<T, -1, 1> layered_3D_2_pillar_1D(Matrix<Matrix<T, -1, -1>, -1, 1> S3)
{
	int nt = (int)S3.size();
	int layer_size = (int)S3(0).size();
	Matrix<T, -1, 1> s((long int)nt*(long int)layer_size);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < nt; i++)
	{
		memcpy(&s(0 + i * layer_size), &(S3(i)(0)), layer_size * sizeof(T));
	}
	return s;
}

//Function for transformations between pillar 2D matrix and pillar 1D vector
template <typename T>
Matrix<T, -1, 1> pillar_2D_2_pillar_1D(Matrix<T, -1, -1> A)
{
	int A_size = (int)A.size();
	Matrix<T, -1, 1> arr(A_size);
	memcpy(&arr(0), &A(0), A_size * sizeof(T));
	return arr;
}
template <typename T>
Matrix<T, -1, -1> pillar_1D_2_pillar_2D(Matrix<T, -1, 1> arr, int row_num, int col_num)
{
	int arr_size = (int)arr.size();
	if (arr_size != row_num * col_num)
	{
		cout << "The size of arr must be equal to row_num*col_num!" << endl;
		exit(1);
	}
	Matrix<T, -1, -1> A(row_num, col_num);
	memcpy(&A(0), &arr(0), arr_size * sizeof(T));
	return A;
}

//Function for getting the maximum amplitude of 'VectorXd' or 'double'
double max_amp(VectorXd x)
{
	return x.cwiseAbs().maxCoeff();
}
double max_amp(double x)
{
	return abs(x);
}

//Function for getting the norm of 'VectorXd' or 'double'
double norm(VectorXd x)
{
	return x.norm();
}
double norm(double x)
{
	return abs(x);
}

//Function for getting the norm of high-dimensional matrix
template <typename T>
double norm(T S)
{
	double S_energy = 0;
	int i;
	for (i = 0; i < S.size(); i++)
	{
		S_energy += S(i).squaredNorm();
	}
	return sqrt(S_energy);
}

//Function for checking the state of solver 
template <typename T>
void check_solver_state(T &solver, int check_type)
{
	string prompt;
	switch (check_type)
	{
	case 1:
		prompt = solver.info() == Success ? "Eigen factorization succeeded." :
			"Eigen factorization failed!";
		break;
	case 2:
		prompt = solver.info() == Success ? "Eigen solving succeeded." : "Eigen solving failed!";
		break;
	default:
		cout << "check_type is out of range!" << endl;
		exit(1);
	}
	cout << prompt << endl;
}

//Function for writting matrix into a text file
template <typename T>
void write_txt(T data, const char *filename)
{
	//Initializaion
	int row_num = (int)data.rows();
	int col_num = (int)data.cols();
	ofstream out;
	out.open(filename);
	if (!out)
	{
		cout << "write_txt cannot open file: " << filename << endl;
		exit(1);
	}

	//Write text file
	int i, j;
	for (j = 0; j < row_num; j++)
	{
		for (i = 0; i < col_num; i++)
		{
			out << data.coeffRef(j, i) << " ";
		}
		out << endl;
	}

	//Close file
	out.close();
	out.clear();
}

//Function for writting matrix into a csv file
template <typename T>
void write_csv(T data, double dt, const char *filename)
{
	//Initialization
	int row_num = (int)data.rows();
	int col_num = (int)data.cols();
	ofstream out;
	out.open(filename);
	if (!out)
	{
		cout << "Cannot open the file!" << endl << filename << endl << endl;
		exit(1);
	}
	int i, j;

	//Write csv file
	for (j = 0; j < row_num; j++)
	{
		out << j * dt << ',';
		if (j < col_num)
		{
			out << j + 1;
		}
		for (i = 0; i < col_num; i++)
		{
			out << ',' << data(j, i);
		}
		out << endl;
	}

	//Close file
	out.close();
	out.clear();
}

//Function for getting the null space of the matrix
MatrixXd get_null_space(MatrixXd A)
{
	int n = (int)A.cols();
	int r = (int)A.jacobiSvd().rank();
	return r == n ? MatrixXd(n, 0) : A.jacobiSvd(ComputeFullV).matrixV().rightCols(n - r);
}

//Function for reading one-line matrix
MatrixXd read_one_line_mat(const char *filename)
{
	//Initialization
	ifstream in;
	in.open(filename);
	if (!in)
	{
		cout << "read_one_line_mat cannot open file: " << filename << endl;
		exit(1);
	}
	char separator;
	int i, j;

	//Get the size of matrix
	int row_num, col_num;
	string str;
	getline(in, str);
	int comma_num = 0, semicolon_num = 0;
	for (i = 0; i < str.size(); i++)
	{
		if (str[i] == ',')
		{
			comma_num++;
		}
		else
		{
			if (str[i] == ';')
			{
				semicolon_num++;
			}
		}
	}
	row_num = semicolon_num + 1;
	col_num = comma_num / row_num + 1;
	str.clear();

	//Read matrix
	in.seekg(1);
	MatrixXd A(row_num, col_num);
	for (j = 0; j < row_num; j++)
	{
		for (i = 0; i < col_num; i++)
		{
			in >> A(j, i) >> separator;
		}
	}


	//Close file
	in.close();
	in.clear();

	return A;
}

//Function for writting one-line matrix
void write_one_line_mat(MatrixXd A, int digit_num, char decimal_point_sign, bool use_effect_digt_num,
	bool del_redundant_zero, const char* filename)
{
	//Initialization
	ofstream out;
	out.open(filename);
	if (!out)
	{
		cout << "write_one_line_mat cannot open file: " << filename << endl;
		exit(1);
	}
	int i, j;

	//Write matrix
	int row_num = (int)A.rows();
	int col_num = (int)A.cols();
	out << '[';
	for (j = 0; j < row_num; j++)
	{
		for (i = 0; i < col_num - 1; i++)
		{
			out << double2str(A(j, i), digit_num, decimal_point_sign, use_effect_digt_num,
				del_redundant_zero) << ", ";
		}
		out << double2str(A(j, i), digit_num, decimal_point_sign, use_effect_digt_num,
			del_redundant_zero);
		if (j < row_num - 1)
		{
			out << "; ";
		}
		else
		{
			out << ']';
		}
	}

	//Close file
	out.close();
	out.clear();
}

//Function for 1D trapezoidal integral
double trapezoid_integral(RowVectorXd s, double grid_size)
{
	s(1) /= 2;
	s(s.size() - 1) /= 2;
	return s.sum()*grid_size;
}

//Function for 2D trapezoidal integral
double trapezoid_integral2(MatrixXd S, double grid_area)
{
	int row_num = (int)S.rows();
	int col_num = (int)S.cols();
	S.row(0) /= 2;
	S.row(row_num - 1) /= 2;
	S.col(0) /= 2;
	S.col(col_num - 1) /= 2;
	return S.sum()*grid_area;
}

//Function for getting upper triangular of symmetric approximation
template <typename T>
SparseMatrix<T, RowMajor> get_sym_appro_upper_triangular(SparseMatrix<T, RowMajor> &A)
{
	RowVectorXi row_ind, col_ind;
	Matrix<T, 1, -1> val;
	sp2coo(A, row_ind, col_ind, val);
	int nnz = (int)A.nonZeros();
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < nnz; i++)
	{
		if (row_ind(i) > col_ind(i))
		{
			exchange(row_ind(i), col_ind(i));
		}
		else
		{
			if (row_ind(i) == col_ind(i))
			{
				val(i) += val(i);
			}
		}
	}
	val /= 2;
	return coo2sp(row_ind, col_ind, val, (int)A.rows(), (int)A.cols());
}

template <typename T>
void gen_Givens_transform(T x1, T x2, T &cos_val, T &sin_val)
{
	if (abs(x2) <= DBL_EPSILON)
	{
		cos_val = 1;
		sin_val = 0;
		return;
	}
	double amp = sqrt(pow_int(abs(x1), 2) + pow_int(abs(x2), 2));
	cos_val = x1 / amp;
	sin_val = x2 / amp;
}

template <typename T>
void apply_Givens_transform(T &x1, T &x2, T &cos_val, T &sin_val)
{
	T tmp = conj(cos_val)*x1 + conj(sin_val)*x2;
	x2 = -sin_val * x1 + cos_val * x2;
	x1 = tmp;
}

template <typename T>
VectorXd get_col_norm(SparseMatrix<T, RowMajor> &A)
{
	A = A.transpose();
	VectorXd col_norm = get_row_norm(A);
	A = A.transpose();
	return col_norm;
}

template <typename T>
VectorXd get_row_norm(SparseMatrix<T, RowMajor> &A)
{
	int *row_ptr = A.outerIndexPtr();
	T *val_ptr = A.valuePtr();
	int row_num = (int)A.rows();
	RowVectorXi j_vec(row_num);
	VectorXd row_norm = VectorXd::Zero(row_num);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < row_num; i++)
	{
		for (j_vec(i) = row_ptr[i]; j_vec(i) < row_ptr[i + 1]; j_vec(i)++)
		{
			row_norm(i) += pow_int(abs(val_ptr[j_vec(i)]), 2);
		}
		row_norm(i) = sqrt(row_norm(i));
	}
	return row_norm;
}

template <typename T>
void apply_row_scaling(SparseMatrix<T, RowMajor> &A, VectorXd row_scaling)
{
	int *row_ptr = A.outerIndexPtr();
	T *val_ptr = A.valuePtr();
	int row_num = (int)A.rows();
	RowVectorXi j_vec(row_num);
	int i;
#pragma omp parallel for private(i)
	for (i = 0; i < row_num; i++)
	{
		for (j_vec(i) = row_ptr[i]; j_vec(i) < row_ptr[i + 1]; j_vec(i)++)
		{
			val_ptr[j_vec(i)] *= row_scaling(i);
		}
	}
}

template <typename T>
void apply_col_scaling(SparseMatrix<T, RowMajor> &A, VectorXd col_scaling)
{
	A = A.transpose();
	apply_row_scaling(A, col_scaling);
	A = A.transpose();
}

template <typename T>
Matrix<T, -1, 1> get_opt_appro(Matrix<T, -1, -1> A, Matrix<T, -1, 1> b)
{
	return A.fullPivHouseholderQr().solve(b);
}

#endif 
