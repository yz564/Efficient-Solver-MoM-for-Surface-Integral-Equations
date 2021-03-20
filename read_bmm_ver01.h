//Functions for reading binary matrix market format file
//////////////////////////////////////////////////////
//filename is the path of the bmm file to be opened
//////////////////////////////////////////////////////
//A is the data matrix to be read
//////////////////////////////////////////////////////
//By Wenhao Xu, Xi'an Jiaotong University & Duke University

#ifndef READ_BMM_VER01_H
#define READ_BMM_VER01_H

#include "Eigen_plus_v1.h"

void read_bmm_info(const char *filename, bool &is_sparse, bool &is_real)
{
	//Initialization
	FILE *fp;
	fp = fopen(filename, "rb");
	if (!fp)
	{
		cout << "\nread_bmm cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}

	//Close file
	fclose(fp);
}

void read_bmm(const char *filename, SparseMatrix<complex<double>, RowMajor> &A)
{
	//Initialization
	FILE *fp;
	fp = fopen(filename, "rb");
	if (!fp)
	{
		cout << "\nread_bmm cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}
	bool is_sparse, is_real;
	fread(&is_sparse, sizeof(bool), 1, fp);
	fread(&is_real, sizeof(bool), 1, fp);
	if (!is_sparse || is_real)
	{
		cout << "read_bmm: the objective bmm file must be in sparse and complex<double> format!" << endl;
		exit(1);
	}
	int row_num, col_num;
	fread(&row_num, sizeof(int), 1, fp);
	fread(&col_num, sizeof(int), 1, fp);
	int nnz;
	fread(&nnz, sizeof(int), 1, fp);

	//Read data
	RowVectorXi row_ind(nnz), col_ind(nnz);
	fread(&row_ind(0), sizeof(int), nnz, fp);
	fread(&col_ind(0), sizeof(int), nnz, fp);
	RowVectorXd real_val(nnz), imag_val(nnz);
	fread(&real_val(0), sizeof(double), nnz, fp);
	fread(&imag_val(0), sizeof(double), nnz, fp);
	RowVectorXcd val(nnz);
	val.real() = real_val;
	val.imag() = imag_val;
	A = coo2sp(row_ind, col_ind, val, row_num, col_num);

	//Close file
	fclose(fp);
}

void read_bmm(const char *filename, MatrixXcd &A)
{
	//Initialization
	FILE *fp;
	fp = fopen(filename, "rb");
	if (!fp)
	{
		cout << "\nread_bmm cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}
	bool is_sparse, is_real;
	fread(&is_sparse, sizeof(bool), 1, fp);
	fread(&is_real, sizeof(bool), 1, fp);
	if (is_sparse || is_real)
	{
		cout << "read_bmm: the objective bmm file must be in dense and complex<double> format!" << endl;
		exit(1);
	}
	int row_num, col_num;
	fread(&row_num, sizeof(int), 1, fp);
	fread(&col_num, sizeof(int), 1, fp);
	long int nnz = (long int)row_num*col_num;

	//Read data
	RowVectorXd real_val(nnz), imag_val(nnz);
	fread(&real_val(0), sizeof(double), nnz, fp);
	fread(&imag_val(0), sizeof(double), nnz, fp);
	RowVectorXcd val(nnz);
	val.real() = real_val;
	val.imag() = imag_val;
	A = MatrixXcd(row_num, col_num);
	memcpy(&A(0), &val(0), nnz * sizeof(complex<double>));

	//Close file
	fclose(fp);
}

void read_bmm(const char *filename, MatrixXd &A)
{
	//Initialization
	FILE *fp;
	fp = fopen(filename, "rb");
	if (!fp)
	{
		cout << "\nread_bmm cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}
	bool is_sparse, is_real;
	fread(&is_sparse, sizeof(bool), 1, fp);
	fread(&is_real, sizeof(bool), 1, fp);
	if (is_sparse || !is_real)
	{
		cout << "read_bmm: the objective bmm file must be in dense and double format!" << endl;
		exit(1);
	}
	int row_num, col_num;
	fread(&row_num, sizeof(int), 1, fp);
	fread(&col_num, sizeof(int), 1, fp);
	long int nnz = (long int)row_num*col_num;

	//Read data
	A = MatrixXd(row_num, col_num);
	fread(&A(0), sizeof(double), nnz, fp);

	//Close file
	fclose(fp);
}

void read_bmm(const char *filename, VectorXi &ordering_vec)
{
	MatrixXd A;
	read_bmm(filename, A);
	MatrixXi ordering_mat = cwise_convert<int>(A);
	int n = (int)ordering_mat.size();
	ordering_vec = VectorXi(n);
	memcpy(&ordering_vec(0), &ordering_mat(0), n * sizeof(int));
}

#endif

