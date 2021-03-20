#ifndef PARAMETERS_H// include guard
#define PARAMETERS_H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "basicmath.h"

std::istream& get_line_strip_comments(std::istream& stm, std::string& str)
{
	if (std::getline(stm, str))
	{
		auto pos = str.find("#");
		if (pos == 0) return get_line_strip_comments(stm, str);
		else if (pos != std::string::npos) str.erase(pos);
	}
	return stm;
}

template<typename T>
void ReadLine_helper(const std::string & str, int num, std::vector<T> & data) {
	std::string word = "";
	int n = 0;
	for (size_t i = 0; i < str.length(); ++i) {
		char x = str[i];
		if (x == '\0' || n >= num) {
			break;
		}
		else if (x != ' ') {
			word += x;
		}
		else if (word.length() > 0) {
			if (std::is_floating_point<T>::value) {
				data[n] = std::stod(word);
			}
			else {
				data[n] = std::stoi(word);
			}
			word = "";
			n++;
		}
	}
	if (n < num && word.length()>0) {
		if (std::is_floating_point<T>::value) {
			data[n] = std::stod(word);
		}
		else {
			data[n] = std::stoi(word);
		}
	}
}

template <typename T>
class InputInfo {
public:
	T freq;
	std::vector<T> eps_r_out;
	std::vector<T> mu_r_out;
	std::vector<T> cond_out;
	std::vector<T> losstangent_out;
	std::vector<T> eps_r_in;
	std::vector<T> mu_r_in;
	std::vector<T> cond_in;
	std::vector<T> losstangent_in;
	int source_num;
	std::vector<int> source_type;
	std::vector<T> source_magnitude;
	std::vector<T> source_phase;
	std::vector<T> source_theta;
	std::vector<T> source_phi;
	std::vector <std::vector<T> > source_location;
	int receiver_numbers;
	T receiver_start[3];
	T receiver_end[3];
	int innerGaussDegree;
	int outerGaussDegree;
	T singularThreshold;
	int allowedthreads;
	std::string mesh_file;

	void print() {
		std::cout << "frequency: " << freq << "Hz\n";
		std::cout << "the number of receivers: " << receiver_numbers << "\n";
		std::cout << "receiver location starts from " << "x=" << receiver_start[0] << " y=" << receiver_start[1] << "z=" << receiver_start[2] << "\n";
		std::cout << "receiver location ends at " << "x=" << receiver_end[0] << " y=" << receiver_end[1] << "z=" << receiver_end[2] << "\n";

		std::ofstream logfile;
		logfile.open("logfile.txt", std::ios_base::app);
		logfile << "frequency: " << freq << "Hz\n";
		logfile << "the number of receivers: " << receiver_numbers << "\n";
		logfile << "receiver location starts from " << "x=" << receiver_start[0] << " y=" << receiver_start[1] << " z=" << receiver_start[2] << "\n";
		logfile << "receiver location ends at " << "x=" << receiver_end[0] << " y=" << receiver_end[1] << " z=" << receiver_end[2] << "\n";
		logfile.close();
	}
};


template <typename T>
InputInfo<T> * ReadInput(std::string inputfile) {
	std::ifstream readfile(inputfile);
	InputInfo<T> * info = new InputInfo<T>();
	std::string str;
	get_line_strip_comments(readfile, str); //ignore the comments start by '#'
	std::vector<T> tmp(5);
	std::vector<int> tmpi(3);
	ReadLine_helper<T>(str, 1, tmp);
	info->freq = tmp[0];
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmpi);
	int MaterialNum = tmpi[0];
	for (int i = 0; i < MaterialNum; ++i) {
		get_line_strip_comments(readfile, str);
		ReadLine_helper<T>(str, 4, tmp);
		(info->eps_r_out).push_back(tmp[0]);
		(info->mu_r_out).push_back(tmp[1]);
		(info->cond_out).push_back(tmp[2]);
		(info->losstangent_out).push_back(tmp[3]);
		get_line_strip_comments(readfile, str);
		ReadLine_helper<T>(str, 4, tmp);
		(info->eps_r_in).push_back(tmp[0]);
		(info->mu_r_in).push_back(tmp[1]);
		(info->cond_in).push_back(tmp[2]);
		(info->losstangent_in).push_back(tmp[3]);
	}
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmpi);
	info->source_num = tmpi[0];
	for (int i = 0; i < info->source_num; ++i) {
		get_line_strip_comments(readfile, str);
		ReadLine_helper<T>(str, 5, tmp);
		(info->source_magnitude).push_back(tmp[0]);
		(info->source_phase).push_back(tmp[1]);
		(info->source_theta).push_back(tmp[2] * PI<T> / 180.0);
		(info->source_phi).push_back(tmp[3] * PI<T> / 180.0);
		(info->source_type).push_back(tmp[4]);
		get_line_strip_comments(readfile, str);
		std::vector<T> tmp_(3);
		ReadLine_helper<T>(str, 3, tmp_);
		(info->source_location).push_back(tmp_);
	}
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmpi);
	info->receiver_numbers = tmpi[0];
	get_line_strip_comments(readfile, str);
	ReadLine_helper<T>(str, 3, tmp);
	info->receiver_start[0] = tmp[0];
	info->receiver_start[1] = tmp[1];
	info->receiver_start[2] = tmp[2];
	get_line_strip_comments(readfile, str);
	ReadLine_helper<T>(str, 3, tmp);
	info->receiver_end[0] = tmp[0];
	info->receiver_end[1] = tmp[1];
	info->receiver_end[2] = tmp[2];

	
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 2, tmpi);
	info->innerGaussDegree = tmpi[0];
	info->outerGaussDegree = tmpi[1];
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmpi);
	info->singularThreshold = tmpi[0];
	get_line_strip_comments(readfile, str);
	ReadLine_helper<int>(str, 1, tmpi);
	info->allowedthreads = tmpi[0];

	get_line_strip_comments(readfile, str);
	info->mesh_file = str;
	return info;
}


#endif