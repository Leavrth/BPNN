#include "DataStruct.h"
using namespace TR;
namespace TR {
	size_t num_e_1 = 0;
	size_t num_e_2 = 24;
	size_t num_e_3 = 12;
	size_t pre_num = 0;
	size_t pre_num_res = 0;
}
fitter::~fitter() {
	delete[] mu;
	delete[] sigma;
}

interval::~interval() {
	delete[] tmin;
	delete[] tmax;
	delete[] pov;
}

data::data(const std::string & line, size_t &cnt, fitter & ft, size_t res_index) {
	e_1 = new double[num_e_1];
	size_t begin = 0, end = -1, index = 0;
	double cot;
	bool nt_slv_res = true;
	for (auto c : line) {
		++end;
		if (c == ',') {
			strTod(line.substr(begin, end - begin), cot);
			if (nt_slv_res && index == res_index) {
				res = log(cot);
				ft.mu_res += (res - ft.mu_res) / cnt;
				nt_slv_res = false;
			}
			else {
				e_1[index] = cot;
				e_1[index + pre_num_res] = log(e_1[index]);
				ft.mu[index] += (e_1[index] - ft.mu[index]) / cnt;
				ft.mu[index + pre_num_res] += (e_1[index + pre_num_res] - ft.mu[index + pre_num_res]) / cnt;
				++index;
			}
			begin = end + 1;
		}
	}
	
	++cnt;
}


data::data(const std::string & line, size_t res_index) {
	e_1 = new double[num_e_1];
	size_t begin = 0, end = -1, index = 0;
	double cot;
	bool nt_slv_res = true;
	for (auto c : line) {
		++end;
		if (c == ',') {
			strTod(line.substr(begin, end - begin), cot);
			if (nt_slv_res && index == res_index) {
				res = log(cot);
				nt_slv_res = false;
			}
			else {
				e_1[index] = cot;
				e_1[index + pre_num_res] = log(e_1[index]);
				++index;
			}
			begin = end + 1;
		}
	}

}
data::data(data && rhs) {
	e_1 = rhs.e_1; rhs.e_1 = nullptr;
	res = rhs.res;
}
data::~data() { delete[] e_1; }

trainres::trainres() {
	W_1 = new double*[num_e_2];
	dW_1 = new double*[num_e_2];
	W_2 = new double*[num_e_3];
	dW_2 = new double*[num_e_3];
	W_3 = new double[num_e_3];
	dW_3 = new double[num_e_3];
	b_1 = new double[num_e_2];
	db_1 = new double[num_e_2];
	b_2 = new double[num_e_3];
	db_2 = new double[num_e_3];
	
	srand((unsigned int)time(0));
	// 初始化
	// Xavier initialization 对 tanh 激活函数权重初始化
	const double Xavier_bound = sqrt(6.0 / (num_e_1 + num_e_2));
	for (int i = 0; i != num_e_2; ++i) {
		W_1[i] = new double[num_e_1];
		dW_1[i] = new double[num_e_1];
		for (int j = 0; j != num_e_1; ++j)
			W_1[i][j] = (rand() / double(RAND_MAX) * 2 * Xavier_bound) - Xavier_bound;
	}
	for (int j = 0; j != num_e_2; ++j)
		b_1[j] = 0;
	// He initialization 对ReLU 激活函数权重初始化
	const double He_bound = sqrt(6.0 / num_e_2);
	for (int i = 0; i != num_e_3; ++i) {
		W_2[i] = new double[num_e_2];
		dW_2[i] = new double[num_e_2];
		for (int j = 0; j != num_e_2; ++j)
			W_2[i][j] = (rand() / double(RAND_MAX) * 2 * He_bound) - He_bound;
	}
	for (int j = 0; j != num_e_3; ++j)
		b_2[j] = 0;
	// liner 初始化
	const double liner_bound = sqrt(6.0 / num_e_3);
	for (int j = 0; j != num_e_3; ++j)
		W_3[j] = (rand() / double(RAND_MAX) * 2 * liner_bound) - liner_bound;
	b_3 = 0;
}

trainres::~trainres() {
	delete[] W_1;
	delete[] W_2;
	delete[] W_3;
	delete[] b_1;
	delete[] b_2;
	delete[] dW_1;
	delete[] dW_2;
	delete[] dW_3;
	delete[] db_1;
	delete[] db_2;
}

void TR::strTod(const std::string &str, double & coe) { std::stringstream sstream; sstream << str; sstream >> coe; }
int TR::indexOf(const std::string &str, char c, size_t i) {
	for (; i < str.length(); ++i)
		if (str.at(i) == c) break;
	return i;
}
double TR::tanh(double z) {
	double res = 0;
	double e_z = 0;
	if (z > 0) {
		e_z = square(exp(-z));
		res = (1 - e_z) / (1 + e_z);
	}
	else {
		e_z = square(exp(z));
		res = (e_z - 1) / (e_z + 1);
	}
	return res;
}
double TR::ReLU(double z) {
	return z > 0 ? z : 0;
}
double TR::dReLU(double z) {
	return z > 0 ? 1 : 0;
}
double TR::square(double z) {
	return z * z;
}
bool TR::isna(const std::string &str) {
	bool ret = str.size() > 0 ? false : true;
	for (auto c : str) {
		if ((c >= '0' && c <= '9') || c == '.') continue;
		ret = true; break;
	}
	return ret;
}
void TR::memset2(double **dst, size_t _1, size_t _2, double val) {
	for (int i = 0; i != _1; ++i)
		for (int j = 0; j != _2; ++j)
			dst[i][j] = val;
}
std::string *TR::build_num() {
	std::string *res = new std::string[pre_num_res];
	for (int i = 0; i != pre_num_res; ++i) {
		std::stringstream sstream;
		sstream << i;
		sstream >> res[i];
	}
	return res;
}
