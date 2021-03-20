//Functions for reading binary matrix market format file
//////////////////////////////////////////////////////
//A is the data matrix to be written
//filename is the path of the bmm file to be opened
//////////////////////////////////////////////////////
//By Wenhao Xu, Xi'an Jiaotong University & Duke University
#ifndef WRITE_BMM_VER01_H
#define WRITE_BMM_VER01_H

#include "Eigen_plus_v1.h"

void write_bmm(SparseMatrix<complex<double>, RowMajor> A, const char *filename)
{

	//Initialization
	FILE *fp;
	fp = fopen(filename, "wb");
	if (!fp)
	{
		cout << "\nwrite_sgy cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}
	bool is_sparse = true;
	bool is_real = false;
	int row_num = (int)A.rows();
	int col_num = (int)A.cols();
	int nnz = (int)A.nonZeros();
	RowVectorXi row_ind, col_ind;
	RowVectorXcd val;
	sp2coo(A, row_ind, col_ind, val);
	RowVectorXd real_val=val.real(), imag_val=val.imag();

	//Write data
	fwrite(&is_sparse, sizeof(bool), 1, fp);
	fwrite(&is_real, sizeof(bool), 1, fp);
	fwrite(&row_num, sizeof(int), 1, fp);
	fwrite(&col_num, sizeof(int), 1, fp);
	fwrite(&nnz, sizeof(int), 1, fp);
	fwrite(&row_ind(0), sizeof(int), nnz, fp);
	fwrite(&col_ind(0), sizeof(int), nnz, fp);
	fwrite(&real_val(0), sizeof(double), nnz, fp);
	fwrite(&imag_val(0), sizeof(double), nnz, fp);

	//Close file
	fclose(fp);
}

void write_bmm(MatrixXcd A, const char *filename)
{

	//Initialization
	FILE *fp;
	fp = fopen(filename, "wb");
	if (!fp)
	{
		cout << "\nwrite_sgy cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}
	bool is_sparse = false;
	bool is_real = false;
	int row_num = (int)A.rows();
	int col_num = (int)A.cols();
	MatrixXd real_A = A.real();
	MatrixXd imag_A = A.imag();
	long int nnz = (long int)row_num*col_num;

	//Write data
	fwrite(&is_sparse, sizeof(bool), 1, fp);
	fwrite(&is_real, sizeof(bool), 1, fp);
	fwrite(&row_num, sizeof(int), 1, fp);
	fwrite(&col_num, sizeof(int), 1, fp);
	fwrite(&real_A(0), sizeof(double), nnz, fp);
	fwrite(&imag_A(0), sizeof(double), nnz, fp);
	
	//Close file
	fclose(fp);
}

void write_bmm(MatrixXd A, const char *filename)
{

	//Initialization
	FILE *fp;
	fp = fopen(filename, "wb");
	if (!fp)
	{
		cout << "\nwrite_sgy cannot open file: " << endl << filename << endl << endl;
		exit(1);
	}
	bool is_sparse = false;
	bool is_real = true;
	int row_num = (int)A.rows();
	int col_num = (int)A.cols();
	long int nnz = (long int)row_num*col_num;

	//Write data
	fwrite(&is_sparse, sizeof(bool), 1, fp);
	fwrite(&is_real, sizeof(bool), 1, fp);
	fwrite(&row_num, sizeof(int), 1, fp);
	fwrite(&col_num, sizeof(int), 1, fp);
	fwrite(&A(0), sizeof(double), nnz, fp);

	//Close file
	fclose(fp);
}

#endif