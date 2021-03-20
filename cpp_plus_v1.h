//Program for adding some basic functions for C++
/////////////////////////////////////////////////////////////
//By Wenhao Xu, Xi'an Jiaotong University
#pragma once
#include <algorithm>
#include <complex>
#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;

//Add class for 2D matrix
template<typename T>
class simple_matrix
{
public:
	int  row_num;
	int  col_num;
	T* p;
public:
	simple_matrix(int row_num, int col_num)
	{
		this->row_num = row_num;
		this->col_num = col_num;
		this->p = new T[row_num*col_num];
		memset(p, 0, row_num*col_num * sizeof(T));
	}
	T& operator ()(int row_ind, int col_ind)
	{
		return *(p + col_ind * row_num + row_ind);
	}
};

//Add class for binary tree
template <typename T>
class simple_binary_tree
{
public:
	T element;
	simple_binary_tree<T>* left_node;
	simple_binary_tree<T>* right_node;

public:
	static simple_binary_tree<T> *preorder_create(T NULL_flag)
	{
		T val;
		cin >> val;
		if (val == NULL_flag)
		{
			return NULL;
		}
		else
		{
			simple_binary_tree<T> *current_node = new simple_binary_tree<T>;
			current_node->element = val;
			current_node->left_node = preorder_create(NULL_flag);
			current_node->right_node = preorder_create(NULL_flag);
			return current_node;
		}
	}
	void preorder_print()
	{
		preorder_print(this);
		cout << endl;
	}
	static void preorder_print(simple_binary_tree<T> * root)
	{
		if (root != NULL)
		{
			cout << root->element << ' ';
			preorder_print(root->left_node);
			preorder_print(root->right_node);
		}
	}
	void inorder_print()
	{
		inorder_print(this);
		cout << endl;
	}
	static void inorder_print(simple_binary_tree<T> * root)
	{
		if (root != NULL)
		{
			inorder_print(root->left_node);
			cout << root->element << ' ';
			inorder_print(root->right_node);
		}
	}
};

//Add function for judging Linux system
bool judge_Linux_system()
{
#if defined(__linux__)
	return true;
#else
	return false;
#endif
}

//Add function for judging Linux system
bool judge_windows_system()
{
#if defined(_WIN32)
	return true;
#else
	return false;
#endif
}

//Add function for printing vector
template <typename T>
void print_vec(T* x, int len)
{

	//Declaration
	int i;

	//Do the printing
	cout << endl << '[';
	for (i = 0; i < len - 1; i++)
	{
		cout << x[i] << " ";
	}
	cout << x[i] << ']' << endl << endl;
}


//Add function for exchanging values
template <typename T>
void exchange(T &a, T &b)
{
	T tmp;
	tmp = a;
	a = b;
	b = tmp;
}

//Get left memory (unit of MB)
double get_Linux_left_mem()
{
	int mem_total = -1, mem_free = -1, mem_buffers = -1, mem_cached = -1;
	char name[20];
	FILE *fp;
	char buf1[128], buf2[128], buf3[128], buf4[128], buf5[128];
	int buff_len = 128;
	fp = fopen("/proc/meminfo", "r");
	if (fp == NULL)
	{
		std::cerr << "get_left_mem() error! File not exist" << std::endl;
		exit(-1);
	}
	if (NULL == fgets(buf1, buff_len, fp) ||
		NULL == fgets(buf2, buff_len, fp) ||
		NULL == fgets(buf3, buff_len, fp) ||
		NULL == fgets(buf4, buff_len, fp) ||
		NULL == fgets(buf5, buff_len, fp))
	{
		std::cerr << "get_left_mem() error! Fail to read!" << std::endl;
		fclose(fp);
		exit(-1);
	}
	fclose(fp);
	sscanf(buf1, "%s%d", name, &mem_total); //Unit of KB
	sscanf(buf2, "%s%d", name, &mem_free);
	sscanf(buf4, "%s%d", name, &mem_buffers);
	sscanf(buf5, "%s%d", name, &mem_cached);
	return double(mem_free + mem_buffers + mem_cached) / 1024;
}

//Add function for getting safe index
int safe_ind(int ind, int max_ind_plus_1, int min_ind = 0)
{
	return max(min(ind, max_ind_plus_1 - 1), min_ind);
}

//Add function for replacing all distinct substrings of a string
string & replace_all_distinct(string & str, const string & old_sub_str,
	const string & new_sub_str)
{
	for (string::size_type pos(0); pos != string::npos; pos += new_sub_str.length())
	{
		if ((pos = str.find(old_sub_str, pos)) != string::npos)
			str.replace(pos, old_sub_str.length(), new_sub_str);
		else   break;
	}
	return str;
}

//Add function for getting maximum adjacent sum
/////////////////////////////////////////////////////////
//x is the vector to produce maxium adjacent sum
//n is the size of x
//m is the length of adjacent elements
/////////////////////////////////////////////////////////
//start_pos is the starting position of returned maximum adjacent sum
//max_sum is the returned maximum adjacent sum
template <typename T>
T get_max_adjacent_sum(T *x, int n, int m, int &start_pos)
{
	//Initialization
	T tmp_sum = x[0];
	int i;
	for (i = 1; i < m; i++)
	{
		tmp_sum += x[i];
	}
	T max_sum = tmp_sum;
	start_pos = 0;

	//Get maximum adjacent sum
	for (i = m; i < n; i++)
	{
		tmp_sum = tmp_sum + x[i] - x[i - m];
		if (tmp_sum > max_sum)
		{
			start_pos = i - m + 1;
			max_sum = tmp_sum;
		}
	}

	return max_sum;
}

//Add function for calculating integral power by multiple multiplications or divisions,
// which should be more accurate than built-in function 'pow' in such case
template <typename T>
T pow_int(T x, int n)
{
	int i;
	T result = 1;
	if (n > 0)
	{
		for (i = 0; i < n; i++)
		{
			result *= x;
		}
	}
	if (n < 0)
	{
		for (i = n; i < 0; i++)
		{
			result /= x;
		}
	}
	return result;
}

//Add functions for counting elements whose amplitudes are greater than or equal to threshold
int count_valid_elem(double *x, int n, double threshold)
{
	int i;
	int count = 0;
	for (i = 0; i < n; i++)
	{
		if (fabs(x[i]) >= threshold)
		{
			count++;
		}
	}
	return count;
}
int count_valid_elem(complex<double> *x, int n, double threshold)
{
	int i;
	int count = 0;
	for (i = 0; i < n; i++)
	{
		if (abs(x[i]) >= threshold)
		{
			count++;
		}
	}
	return count;
}

//Add functions for clearing elements whose amplitudes are greater than threshold
void clear_invalid_elem(double *x, int n, double threshold)
{
	int i;
	for (i = 0; i < n; i++)
	{
		if (fabs(x[i]) <= threshold)
		{
			x[i] = 0;
		}
	}
}
void clear_invalid_elem(complex<double> *x, int n, double threshold)
{
	int i;
	for (i = 0; i < n; i++)
	{
		if (abs(x[i]) <= threshold)
		{
			x[i] = 0;
		}
	}
}

//Add sign function
template <typename T>
T sign(T x)
{
	if (x > 0)
	{
		return 1;
	}
	else if (x < 0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

//Add functions for reversing array 
template<typename T>
void reverse(T *in, T *out, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		out[i] = in[n - 1 - i];
	}
}
template <typename T>
void reverse_in_place(T *x, int n)
{
	T *x1 = new T[n];
	int i;
	for (i = 0; i < n; i++)
	{
		x1[i] = x[n - 1 - i];
	}
	memcpy(x, x1, n * sizeof(T));
	delete x1;
}

//Add factorial function
double factorial(int n)
{
	double result = 1;
	for (int i = 2; i <= n; i++)
	{
		result *= i;
	}
	return result;
}

//Add funtion for sorting by key
template <typename T1, typename T2>
class key_val_class
{
public:
	T1 key;
	T2 val;
public:
	bool operator <(const key_val_class<T1, T2> &b) const
	{
		return key < b.key;
	}
	bool operator >(const key_val_class<T1, T2> &b) const
	{
		return key > b.key;
	}
};
template<typename T1, typename T2>
void sort_by_key(T1 *key, T2 *val, int n, bool is_ascend = true)
{
	key_val_class<T1, T2> *A = new key_val_class<T1, T2>[n];
	int i;
	for (i = 0; i < n; i++)
	{
		A[i].key = *(key + i);
		A[i].val = *(val + i);
	}
	if (is_ascend)
	{
		sort(&(A[0]), &(A[0]) + n, less<key_val_class<T1, T2> >());
	}
	else
	{
		sort(&(A[0]), &(A[0]) + n, greater<key_val_class<T1, T2> >());
	}
	for (i = 0; i < n; i++)
	{
		*(key + i) = A[i].key;
		*(val + i) = A[i].val;
	}
	delete A;
}

//Add function for gettting supremum index of specific element by dichotomy
template <typename T>
int get_supremum_ind_by_dichotomy(T *x, int n, T elem)
{
	//Check input
	if (x[0] > x[1] || x[n - 1] < x[n - 2])
	{
		cout << "get_supremum_ind_by_dichotomy:" << endl;
		cout << "x must be sorted in ascending order!" << endl;
		exit(1);
	}

	//Use dichotomy
	int left_ind = 0, right_ind = n - 1;
	int middle_ind;
	while (right_ind - left_ind > 1)
	{
		middle_ind = (left_ind + right_ind) / 2;
		if (T[middle_ind] >= elem)
		{
			right_ind = middle_ind;
		}
		else
		{
			left_ind = middle_ind;
		}
	}

	return right_ind;
}

//Add function for getting the argument of a complex
//(The range of argument is [-180 180])
void get_complex_argument(complex<double> *c, double *argument, int n)
{
	int i;
	double pi = acos(-1.0);
	for (i = 0; i < n; i++)
	{
		if (abs(c[i]) <= DBL_EPSILON)
		{
			argument[i] = 0;
		}
		else
		{
			argument[i] = acos(real(c[i]) / abs(c[i])) / pi * 180;
			if (imag(c[i]) < 0)
			{
				argument[i] = -argument[i];
			}
		}
	}
}

//Add Ricker function
double Ricker_fun(double f0, double t0, double t)
{
	double pi = acos(-1.0);
	return (1 - 2 * pow_int(pi*f0*(t - t0), 2))*exp(-pow_int(pi*f0*(t - t0), 2));
}

//Add function for writting binary file
template <typename T>
void write_bin(const char* filename, T *x, int n)
{
	//Initialization
	FILE *fp;
	fp = fopen(filename, "wb+");
	if (!fp)
	{
		cout << "write_bin cannot open object file!" << endl;
		return;
	}

	//Write binary file
	fwrite(x, sizeof(T), n, fp);

	//Close file
	fclose(fp);
}

//Add function for juedging prime number
template <typename T>
bool is_prime(T n)
{
	if (n == 2 || n == 3)
	{
		return true;
	}
	if (n <= 1 || n % 2 == 0 || n % 3 == 0)
	{
		return false;
	}
	T N = (T)sqrt(n);
	T i;
	for (i = 5; i <= N; i++)
	{
		if (n%i == 0)
		{
			return false;
		}
	}
	return true;
}

//Add function for getting the binary size of a file
int get_file_size(FILE* fp)
{
	int current_pos, len;
	current_pos = ftell(fp);
	fseek(fp, 0, SEEK_END);
	len = ftell(fp);
	fseek(fp, current_pos, SEEK_SET);
	return len;
}

//Add function for switching endian between big endian and little endian
template <typename T>
void switch_endian(T *x, int n)
{
	int i, j;
	T y;
	char *ptr_x, *ptr_y = (char *)&y;
	int type_size = sizeof(T);
	for (i = 0; i < n; i++)
	{
		ptr_x = (char *)(x + i);
		for (j = 0; j < type_size; j++)
		{
			*(ptr_y + j) = *(ptr_x + type_size - 1 - j);
		}
		*(x + i) = y;
	}
}

//Add function for converting IBM float memory into IEEE float memory, based on Fu Maosong (2011)
unsigned long IBM2IEEE(unsigned long float_IBM)
{
	if ((float_IBM << 1) == 0)
	{
		return float_IBM;
	}

	//Initialization
	unsigned long MASK_S = 0x80000000;
	unsigned long MASK_IBM32_E = 0x7F000000;
	unsigned long MASK_IBM32_F = 0x00FFFFFF;

	//Get the sign part of IEEE float
	unsigned long S_IBM_32 = float_IBM & MASK_S;

	//Get the exponent part of IEEE float
	unsigned long E_IBM_32 = float_IBM & MASK_IBM32_E;

	//Get the fractional part of IEEE float
	unsigned long F_IBM_32 = float_IBM & MASK_IBM32_F;

	//Get the fractinal part of IEEE float
	unsigned long radix = 0;
	unsigned long F_IEEE_32 = F_IBM_32;
	while (radix <= 3 && F_IEEE_32 < 0x01000000)
	{
		radix++;
		F_IEEE_32 = F_IEEE_32 << 1;
	}
	F_IEEE_32 = (F_IEEE_32 - 0x01000000) >> 1;

	//Get the exponent part of IEEE float
	unsigned long E_IEEE_32 = (((E_IBM_32 >> 22) - 130) - (radix - 1)) << 23;

	//Process overflow
	if (E_IEEE_32 > 0x7F800000)
	{
		return S_IBM_32 | 0x7F800000;
	}
	if (E_IEEE_32 < 0x10000000)
	{
		return S_IBM_32;
	}
	else
	{
		return S_IBM_32 | E_IEEE_32 | F_IEEE_32;
	}
}

//Add function for converting IEEE float memory into IBM float memory, based on Fu Maosong (2011)
unsigned long IEEE2IBM(unsigned long float_IEEE)
{
	if ((float_IEEE << 1) == 0)
	{
		return float_IEEE;
	}

	//Initialization
	unsigned long MASK_S = 0x80000000;
	unsigned long MASK_IEEE32_E = 0x7F800000;
	unsigned long MASK_IEEE32_F = 0x007FFFFF;
	unsigned long MASK_IBM32_E = 0x7F000000;

	//Get the sign part of IEEE float
	unsigned long S_IEEE_32 = float_IEEE & MASK_S;

	//Get the exponent part of IEEE float
	unsigned long E_IEEE_32 = float_IEEE & MASK_IEEE32_E;

	//Get the fractional part of IEEE float
	unsigned long F_IEEE_32 = float_IEEE & MASK_IEEE32_F;

	//Get the exponent part of IBM float
	unsigned long E_IBM_32 = ((E_IEEE_32 + 0x41000000) >> 1) &MASK_IBM32_E;

	//Get the fractional part of IBM float
	unsigned F_IBM_32 = ((F_IEEE_32 << 1) + 0x1000000) >> 1;

	//When (E_IEEE_32+130)%4!=0, F_IBM_32 should be divided by (E_IEEE_32+130)%4
	unsigned long reminder = (E_IEEE_32 + 0x41000000) >> 23;
	if ((reminder & 3) != 0)
	{
		E_IBM_32 = E_IBM_32 + 0x01000000;
		F_IBM_32 = F_IBM_32 >> (4 - (reminder & 3));
	}

	return S_IEEE_32 | E_IBM_32 | F_IBM_32;
}

//Function for conversion between IEEE float and IBM float, based on Fu Maosong (2011)
//(Conversion_type is the type of conversion, 1 means converting IBM float into IEEE float, 
// 2 means converting IEEE float into IBM float)
void float_conversion(float *buffer, int n, int conversion_type)
{
	//Initialization
	int i;
	unsigned long * float_p = NULL;
	unsigned long float_tmp;

	switch (conversion_type)
	{
	case 1:
		for (i = 0; i < n; i++)
		{
			float_p = (unsigned long*)(buffer + i);
			float_tmp = IBM2IEEE(*float_p);
			buffer[i] = *((float*)(&float_tmp));
		}
		break;
	case 2:
		for (i = 0; i < n; i++)
		{
			float_p = (unsigned long*)(buffer + i);
			float_tmp = IEEE2IBM(*float_p);
			buffer[i] = *((float*)(&float_tmp));
		}
		break;
	default:
		cout << "The conversion_type is out of range!" << endl;
		exit(1);
		break;
	}
}

//Add function for getting the minimum element of three
template <typename T>
T min3(T x, T y, T z)
{
	return min(min(x, y), z);
}

//Add function for getting the maximum element of three
template <typename T>
T max3(T x, T y, T z)
{
	return max(max(x, y), z);
}

//Add function for safe fmod that reduces the error caused by datatype inaccuracy (like fmod(1.2,0.1))
template <typename T>
T safe_fmod(T x, T y)
{
	T r = fmod(x, y);
	return abs(abs(r) - abs(y)) <= DBL_EPSILON ? 0 : r;
}

template <typename T1, typename T2>
T1 simple_smooth(T1 x, void smooth_fun(T2 *, T2 *, int))
{
	int n = (int)x.size();
	T1 y = x;
	smooth_fun(&(x(0)), &(y(0)), n);
	return y;
}

template <typename T>
void lin_smooth3(T *in, T *out, int N)
{
	int i;
	if (N < 3)
	{
		for (i = 0; i <= N - 1; i++)
		{
			out[i] = in[i];
		}
	}
	else
	{
		out[0] = (5.0 * in[0] + 2.0 * in[1] - in[2]) / 6.0;

		for (i = 1; i <= N - 2; i++)
		{
			out[i] = (in[i - 1] + in[i] + in[i + 1]) / 3.0;
		}

		out[N - 1] = (5.0 * in[N - 1] + 2.0 * in[N - 2] - in[N - 3]) / 6.0;
	}
}

template <typename T>
void lin_smooth5(T *in, T *out, int N)
{
	int i;
	if (N < 5)
	{
		for (i = 0; i <= N - 1; i++)
		{
			out[i] = in[i];
		}
	}
	else
	{
		out[0] = (3.0 * in[0] + 2.0 * in[1] + in[2] - in[4]) / 5.0;
		out[1] = (4.0 * in[0] + 3.0 * in[1] + 2.0 * in[2] + in[3]) / 10.0;
		for (i = 2; i <= N - 3; i++)
		{
			out[i] = (in[i - 2] + in[i - 1] + in[i] + in[i + 1] + in[i + 2]) / 5.0;
		}
		out[N - 2] = (4.0 * in[N - 1] + 3.0 * in[N - 2] + 2.0 * in[N - 3] + in[N - 4]) / 10.0;
		out[N - 1] = (3.0 * in[N - 1] + 2.0 * in[N - 2] + in[N - 3] - in[N - 5]) / 5.0;
	}
}

template <typename T>
void lin_smooth7(T *in, T *out, int N)
{
	int i;
	if (N < 7)
	{
		for (i = 0; i <= N - 1; i++)
		{
			out[i] = in[i];
		}
	}
	else
	{
		out[0] = (13.0 * in[0] + 10.0 * in[1] + 7.0 * in[2] + 4.0 * in[3] +
			in[4] - 2.0 * in[5] - 5.0 * in[6]) / 28.0;

		out[1] = (5.0 * in[0] + 4.0 * in[1] + 3.0 * in[2] + 2.0 * in[3] +
			in[4] - in[6]) / 14.0;

		out[2] = (7.0 * in[0] + 6.0 * in[1] + 5.0 * in[2] + 4.0 * in[3] +
			3.0 * in[4] + 2.0 * in[5] + in[6]) / 28.0;

		for (i = 3; i <= N - 4; i++)
		{
			out[i] = (in[i - 3] + in[i - 2] + in[i - 1] + in[i] + in[i + 1] + in[i + 2] + in[i + 3]) / 7.0;
		}

		out[N - 3] = (7.0 * in[N - 1] + 6.0 * in[N - 2] + 5.0 * in[N - 3] +
			4.0 * in[N - 4] + 3.0 * in[N - 5] + 2.0 * in[N - 6] + in[N - 7]) / 28.0;

		out[N - 2] = (5.0 * in[N - 1] + 4.0 * in[N - 2] + 3.0 * in[N - 3] +
			2.0 * in[N - 4] + in[N - 5] - in[N - 7]) / 14.0;

		out[N - 1] = (13.0 * in[N - 1] + 10.0 * in[N - 2] + 7.0 * in[N - 3] +
			4.0 * in[N - 4] + in[N - 5] - 2.0 * in[N - 6] - 5.0 * in[N - 7]) / 28.0;
	}
}

template <typename T>
void quad_smooth5(T *in, T *out, int N)
{
	int i;
	if (N < 5)
	{
		for (i = 0; i <= N - 1; i++)
		{
			out[i] = in[i];
		}
	}
	else
	{
		out[0] = (31.0 * in[0] + 9.0 * in[1] - 3.0 * in[2] - 5.0 * in[3] + 3.0 * in[4]) / 35.0;
		out[1] = (9.0 * in[0] + 13.0 * in[1] + 12.0 * in[2] + 6.0 * in[3] - 5.0 *in[4]) / 35.0;
		for (i = 2; i <= N - 3; i++)
		{
			out[i] = (-3.0 * (in[i - 2] + in[i + 2]) +
				12.0 * (in[i - 1] + in[i + 1]) + 17.0 * in[i]) / 35.0;
		}
		out[N - 2] = (9.0 * in[N - 1] + 13.0 * in[N - 2] + 12.0 * in[N - 3] + 6.0 * in[N - 4] - 5.0 * in[N - 5]) / 35.0;
		out[N - 1] = (31.0 * in[N - 1] + 9.0 * in[N - 2] - 3.0 * in[N - 3] - 5.0 * in[N - 4] + 3.0 * in[N - 5]) / 35.0;
	}
}

template <typename T>
void quad_smooth7(T *in, T *out, int N)
{
	int i;
	if (N < 7)
	{
		for (i = 0; i <= N - 1; i++)
		{
			out[i] = in[i];
		}
	}
	else
	{
		out[0] = (32.0 * in[0] + 15.0 * in[1] + 3.0 * in[2] - 4.0 * in[3] -
			6.0 * in[4] - 3.0 * in[5] + 5.0 * in[6]) / 42.0;

		out[1] = (5.0 * in[0] + 4.0 * in[1] + 3.0 * in[2] + 2.0 * in[3] +
			in[4] - in[6]) / 14.0;

		out[2] = (1.0 * in[0] + 3.0 * in[1] + 4.0 * in[2] + 4.0 * in[3] +
			3.0 * in[4] + 1.0 * in[5] - 2.0 * in[6]) / 14.0;
		for (i = 3; i <= N - 4; i++)
		{
			out[i] = (-2.0 * (in[i - 3] + in[i + 3]) +
				3.0 * (in[i - 2] + in[i + 2]) +
				6.0 * (in[i - 1] + in[i + 1]) + 7.0 * in[i]) / 21.0;
		}
		out[N - 3] = (1.0 * in[N - 1] + 3.0 * in[N - 2] + 4.0 * in[N - 3] +
			4.0 * in[N - 4] + 3.0 * in[N - 5] + 1.0 * in[N - 6] - 2.0 * in[N - 7]) / 14.0;

		out[N - 2] = (5.0 * in[N - 1] + 4.0 * in[N - 2] + 3.0 * in[N - 3] +
			2.0 * in[N - 4] + in[N - 5] - in[N - 7]) / 14.0;

		out[N - 1] = (32.0 * in[N - 1] + 15.0 * in[N - 2] + 3.0 * in[N - 3] -
			4.0 * in[N - 4] - 6.0 * in[N - 5] - 3.0 * in[N - 6] + 5.0 * in[N - 7]) / 42.0;
	}
}

template <typename T>
void cubic_smooth5(T *in, T *out, int N)
{

	int i;
	if (N < 5)
	{
		for (i = 0; i <= N - 1; i++)
			out[i] = in[i];
	}

	else
	{
		out[0] = (69.0 * in[0] + 4.0 * in[1] - 6.0 * in[2] + 4.0 * in[3] - in[4]) / 70.0;
		out[1] = (2.0 * in[0] + 27.0 * in[1] + 12.0 * in[2] - 8.0 * in[3] + 2.0 * in[4]) / 35.0;
		for (i = 2; i <= N - 3; i++)
		{
			out[i] = (-3.0 * (in[i - 2] + in[i + 2]) + 12.0 * (in[i - 1] + in[i + 1]) + 17.0 * in[i]) / 35.0;
		}
		out[N - 2] = (2.0 * in[N - 5] - 8.0 * in[N - 4] + 12.0 * in[N - 3] + 27.0 * in[N - 2] + 2.0 * in[N - 1]) / 35.0;
		out[N - 1] = (-in[N - 5] + 4.0 * in[N - 4] - 6.0 * in[N - 3] + 4.0 * in[N - 2] + 69.0 * in[N - 1]) / 70.0;
	}
	return;
}

template <typename T>
void cubic_smooth7(T *in, T *out, int N)
{
	int i;
	if (N < 7)
	{
		for (i = 0; i <= N - 1; i++)
		{
			out[i] = in[i];
		}
	}
	else
	{
		out[0] = (39.0 * in[0] + 8.0 * in[1] - 4.0 * in[2] - 4.0 * in[3] +
			1.0 * in[4] + 4.0 * in[5] - 2.0 * in[6]) / 42.0;
		out[1] = (8.0 * in[0] + 19.0 * in[1] + 16.0 * in[2] + 6.0 * in[3] -
			4.0 * in[4] - 7.0* in[5] + 4.0 * in[6]) / 42.0;
		out[2] = (-4.0 * in[0] + 16.0 * in[1] + 19.0 * in[2] + 12.0 * in[3] +
			2.0 * in[4] - 4.0 * in[5] + 1.0 * in[6]) / 42.0;
		for (i = 3; i <= N - 4; i++)
		{
			out[i] = (-2.0 * (in[i - 3] + in[i + 3]) +
				3.0 * (in[i - 2] + in[i + 2]) +
				6.0 * (in[i - 1] + in[i + 1]) + 7.0 * in[i]) / 21.0;
		}
		out[N - 3] = (-4.0 * in[N - 1] + 16.0 * in[N - 2] + 19.0 * in[N - 3] +
			12.0 * in[N - 4] + 2.0 * in[N - 5] - 4.0 * in[N - 6] + 1.0 * in[N - 7]) / 42.0;
		out[N - 2] = (8.0 * in[N - 1] + 19.0 * in[N - 2] + 16.0 * in[N - 3] +
			6.0 * in[N - 4] - 4.0 * in[N - 5] - 7.0 * in[N - 6] + 4.0 * in[N - 7]) / 42.0;
		out[N - 1] = (39.0 * in[N - 1] + 8.0 * in[N - 2] - 4.0 * in[N - 3] -
			4.0 * in[N - 4] + 1.0 * in[N - 5] + 4.0 * in[N - 6] - 2.0 * in[N - 7]) / 42.0;
	}
}

//Function for converting an integer to a string
string int2str(int num)
{
	char *p = new char[1000];
	sprintf(p, "%d", num);
	string str = p;
	delete p;
	return str;
}

//Function for converting a string to an integer (negative sign is allowed)
int str2int(const char *str)
{
	return atoi(str);
}

//Function for converting a double to string
//use_effect_digit_num is a bool variable, which decides whether to use effective digit nubmer
string double2str(double x = 3.1415926, int digit_num = 0, char decimal_point_sign = '.',
	bool use_effect_digit_num = false, bool del_redundant_zero = true)
{
	int decimal_point_pos;
	int sign;
	int i;
	string str = use_effect_digit_num ? ecvt(x, digit_num, &decimal_point_pos, &sign) :
		fcvt(x, digit_num, &decimal_point_pos, &sign);
	if (decimal_point_pos <= 0)
	{
		str.insert(0, 1 - decimal_point_pos, '0');
		decimal_point_pos = 1;
	}
	if (decimal_point_pos != str.size())
	{
		str.insert(decimal_point_pos, 1, decimal_point_sign);
	}
	if (sign != 0)
	{
		str.insert(0, 1, '-');
		decimal_point_pos++;
	}
	if (del_redundant_zero) //Delete redundant zero if requiared
	{
		for (i = (int)str.size(); i >= 1; i--)
		{
			if (str[i - 1] != '0')
			{
				break;
			}
		}
		if (i >= decimal_point_pos && i != str.size())
		{
			str.erase(i, str.size() - i);
		}
	}
	if (str[str.size() - 1] == decimal_point_sign) //Delete redundant decimal point sign
	{
		str.erase(str.size() - 1);
	}
	return str;
}

double str2double(const char *str)
{
	return atof(str);
}