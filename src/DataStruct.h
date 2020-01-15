#pragma once
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace TR {
	extern size_t num_e_1;
	extern size_t num_e_2;
	extern size_t num_e_3;
	extern size_t pre_num;
	extern size_t pre_num_res;
#define AXIS_DIV 50
#define PRE_NUM 12
#define PRE_NUM_RES 11
#define NUM_E_1 22
#define NUM_E_2 24
#define NUM_E_3 12

	// zero-mean normalization
	struct fitter {
		double *mu, mu_res;  //mu[NUM_E_1], mu_res;   // 样本均值
		double *sigma, sigma_res;  //sigma[NUM_E_1], sigma_res;   // 样本方差
		~fitter();
	};
	struct interval {
		double *tmin, *tmax; //tmin[PRE_NUM], tmax[PRE_NUM];
		double *pov; //pov[PRE_NUM_RES*(PRE_NUM_RES - 1) / 2];
		/*
			0<=i<=n-1
			(i, j)  ->   [i(n-1) - (i-1)i/2 + j - i - 1]
			(0,1), (0,2), ..., (0, n-1)
			(1,2), (1,3), ..., (1, n-1)
		*/
		~interval();
	};

	struct data {
		double *e_1;//e_1[NUM_E_1];
		double res;  // diameter
		data(const std::string & line, size_t res_index);
		data(const std::string & line, size_t &cnt, fitter & ft, size_t res_index);
		data(data && rhs);
		~data();
	};
	struct trainres {
		double **W_1; //W_1[NUM_E_2][NUM_E_1];
		double **W_2; //W_2[NUM_E_3][NUM_E_2];
		double *W_3; //W_3[NUM_E_3];
		double *b_1, *b_2, b_3; //b_1[NUM_E_2], b_2[NUM_E_3], b_3;
		double **dW_1; //dW_1[NUM_E_2][NUM_E_1];
		double **dW_2; //dW_2[NUM_E_3][NUM_E_2];
		double *dW_3; //dW_3[NUM_E_3];
		double *db_1, *db_2, db_3; //db_1[NUM_E_2], db_2[NUM_E_3], db_3;
		trainres();
		~trainres();
	};
	
	struct Axisdata {
		const char * title;
		const char * X_title, *Y_title;
		double x_min, x_max;
		double y_min, y_max;
		int X_div_num;
		int Y_div_num;
		double *Y;
		Axisdata(const char * title, const char * Xtitle, const char * Ytitle,
			double x_min, double x_max, double y_min, double y_max, int xdn, int ydn, double *y)
			: title(title), X_title(Xtitle), Y_title(Ytitle), x_min(x_min), x_max(x_max),
			y_min(y_min), y_max(y_max), X_div_num(xdn), Y_div_num(ydn), Y(y) {}
		
	};
	
	void strTod(const std::string &str, double & coe);
	int indexOf(const std::string &str, char c, size_t i = 0);
	double tanh(double z);
	double ReLU(double z);
	double dReLU(double z);
	double square(double z);
	bool isna(const std::string &str);
	void memset2(double **dst, size_t _1, size_t _2, double val = 0.0);
	std::string *build_num();
}