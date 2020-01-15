#include "Trainer.h"
using namespace TR;

FileReader::~FileReader() {
	delete[] isnan_cnt;
}

FileReader::FileReader(const std::string & str, const std::string & res_str, size_t line_num, size_t np) : FILE_RD(str.substr(0,str.size()-4)), num_of_partitions(np), res_str(res_str) {
	
	std::ifstream in(str);
	std::string line;
	std::getline(in, line);
	line.append(",");
	size_t begin = 0, end = -1, index = 0;
	size_t res_index = -1; // 结果的索引
	// 读入表头
	for (auto c : line) {
		++end;
		if (c == ',') {
			std::string name = line.substr(begin, end - begin);
			if (res_index == -1 && name == res_str) res_index = index;
			name_list[std::move(name)] = index++;
			begin = end + 1;
		}
	}
	// 读入数据, 避免数据过大，批量处理
	const size_t batch_size = 50000;
	size_t cnt = 0;
	while (cnt < line_num) {
		size_t end_num = cnt + batch_size;
		end_num = end_num > line_num ? line_num : end_num;
		size_t wid = 0;
		//pre_datas.clear();
		std::stringstream sstream; sstream << FILE_RD << "__" << cnt / batch_size << ".csv";
		std::string outfile; sstream >> outfile;
		std::ofstream out(outfile);
		// 筛选
		isnan_cnt = new size_t[name_list.size()]();
		while (cnt++ != end_num && std::getline(in, line)) {
			line.append(",");
			begin = index = 0;
			end = -1;
			//pre_datas.push_back(std::map<size_t, std::string>());
			auto m = std::map<size_t, std::string>();
			bool is_del = false;
			for (auto c : line) {
				++end;
				if (c == ',') {
					if (index == res_index && end == begin) { // 去掉无res的数据
						//pre_datas.pop_back();
						is_del = true;
						break;
					}
					m[index++] = line.substr(begin, end - begin);
					begin = end + 1;
				}
			}
			if (!is_del) {
				out << line << std::endl; ++wid;
				for (size_t i = 0; i != name_list.size(); ++i) {
					if (isna(m[i]) || m[i] == "0") ++isnan_cnt[i];
				}
			}
			

		}

		num_e_1 = 0;
		for (auto &m : name_list) {
			std::cout << m.first << " " << static_cast<double>(isnan_cnt[m.second]) / wid << std::endl;
			if (isnan_cnt[m.second] < isnan_cnt[res_index]) ++num_e_1;
		}

		del(res_index, cnt / batch_size);
		delete[] isnan_cnt;

		out.close();
	}

	pre_num_res = num_e_1;
	pre_num = pre_num_res + 1;
	num_e_1 *= 2;
	in.close();
}

void FileReader::del(size_t res_index, size_t iter_cur_num) {
	std::ofstream out;
	std::stringstream sstream; sstream << FILE_RD << "_" << iter_cur_num << ".csv";
	std::string filename; sstream >> filename;
	out.open(filename.c_str());
	sstream.clear();
	sstream << FILE_RD << "__" << iter_cur_num << ".csv";
	std::string infile; sstream >> infile;
	std::ifstream in(infile);
	FILE_NAMES.emplace(std::move(filename));
	for (auto &p : name_list) {
		if (isnan_cnt[p.second] <= isnan_cnt[res_index])
			out << p.first << ',';
	}
	out << std::endl;
	std::string line;
	while(std::getline(in, line)){
		size_t begin = 0, index = 0, end = -1;
		auto m = std::map<size_t, std::string>();
		for (auto c : line) {
			++end;
			if (c == ',') {
				m[index++] = line.substr(begin, end - begin);
				begin = end + 1;
			}
		}
		bool m_is_nan = false;
		for (auto &p : m) {
			if (isnan_cnt[p.first] <= isnan_cnt[res_index] && isna(p.second)) {
				m_is_nan = true;
				break;
			}
		}
		if (m_is_nan) continue;
		for (auto &p : name_list) {
			if (isnan_cnt[p.second] <= isnan_cnt[res_index])
				out << m[p.second] << ',';
		}
		out << std::endl;
	}
	out.close();
}

int FileReader::hashcode(const std::string &str) {
	unsigned int h = 0;
	char clast = 0;
	for(auto c : str) {
		if (c == ',' && clast == ',') return 0;
			h = c + (h << 5) - h;
			clast = c;
		}
	return h % num_of_partitions + 1; // 0是非训练集，
}

void FileReader::read(std::vector<data> & datas1, std::vector<data> & datas2, fitter & ft, size_t &res_index, std::map<std::string, size_t> &name_list) {
	ft.mu = new double[num_e_1];
	ft.sigma = new double[num_e_1];
	memset(ft.mu, 0, sizeof(double) * num_e_1);
	memset(ft.sigma, 0, sizeof(double) * num_e_1);
	ft.mu_res = ft.sigma_res = 0;
	size_t cnt = 1;
	res_index = -1; // 结果的索引

	for (auto & s : FILE_NAMES) {
		std::ifstream in;
		in.open(s);
		srand((unsigned int)time(0));
		std::string line;
		std::getline(in, line);
		if (res_index == -1) {
			size_t begin = 0, end = -1, index = 0;
			for (auto c : line) {
				++end;
				if (c == ',') {
					std::string name = line.substr(begin, end - begin);
					if (res_index == -1 && name == res_str) res_index = index;
					name_list[std::move(name)] = index++;
					begin = end + 1;
				}
			}
		}
		// loading and caculate mu
		while (std::getline(in, line)) {
			if (rand() / double(RAND_MAX) < 0.8)
				datas1.emplace_back(data(line, cnt, ft, res_index));
			else datas2.emplace_back(data(line, res_index));
		}
	}
	// caculate sigma
	cnt = 1;
	for (auto &v : datas1) {
		for (int i = 0; i != num_e_1; ++i)
			ft.sigma[i] += square(v.e_1[i] - ft.mu[i]) / cnt;
		ft.sigma_res += square(v.res - ft.mu_res) / cnt;
		++cnt;
	}
	for (int i = 0; i != num_e_1; ++i)
		ft.sigma[i] = sqrt(ft.sigma[i] / (cnt - 2) * (cnt - 1));
	ft.sigma_res = sqrt(ft.sigma_res / (cnt - 2) * (cnt - 1));
	
}

void DataTrainer::traindatas(int epoch, int batch_size, double alpha, double lambda) {
	// use Adam optimizer
	const double beta_1 = 0.9;
	double **m_1 = new double*[num_e_2];
	for (int i = 0; i != num_e_2; ++i) m_1[i] = new double[num_e_1]();
	double **m_2 = new double*[num_e_3];
	for (int i = 0; i != num_e_3; ++i) m_2[i] = new double[num_e_2]();
	double *m_3 = new double[num_e_3]();
	double *m_b_1 = new double[num_e_2]();
	double *m_b_2 = new double[num_e_3]();
	double m_b_3 = 0;
	// m_1[NUM_E_2][NUM_E_1], m_2[NUM_E_3][NUM_E_2], m_3[NUM_E_3];
	// m_b_1[NUM_E_2], m_b_2[NUM_E_3], m_b_3;

	// m = beta_1 * m + (1 - beta_1) * delta
	const double beta_2 = 0.999;
	double **v_1 = new double*[num_e_2];
	for (int i = 0; i != num_e_2; ++i) v_1[i] = new double[num_e_1]();
	double **v_2 = new double*[num_e_3];
	for (int i = 0; i != num_e_3; ++i) v_2[i] = new double[num_e_2]();
	double *v_3 = new double[num_e_3]();
	double *v_b_1 = new double[num_e_2]();
	double *v_b_2 = new double[num_e_3]();
	double v_b_3 = 0;
	// v_1[NUM_E_2][NUM_E_1], v_2[NUM_E_3][NUM_E_2], v_3[NUM_E_3];
	// v_b_1[NUM_E_2], v_b_2[NUM_E_3], v_b_3;
	
	// v = beta_2 * v + (1 - beta_2) * delta * delta
	// 避免初期影响，偏差修正，t为递推次数即m_t
	double beta_1_t = beta_1;
	double beta_2_t = beta_2;
	// m_hat = m / (1 - beta_1 ^ t)
	// v_hat = v / (1 - beta_2 ^ t)
	const double epsilon = 1e-8; // 避免除数为0
	// 每次迭代 theta = theta - alpha * m_hat / (sqrt(v_hat) + epsilon)

	double *E_2 = new double[num_e_2];
	double *E_3 = new double[num_e_3];
	double *delta_1 = new double[num_e_2];
	double *delta_2 = new double[num_e_3];
	double delta_3;
	//E_2[NUM_E_2], E_3[NUM_E_3];
	//delta_1[NUM_E_2], delta_2[NUM_E_3], delta_3;
	double y_hat;
	for (int epc = 0; epc != epoch; ++epc) {
		for (unsigned int num_batch = 0; num_batch < train_datas.size(); num_batch += batch_size) {
			unsigned int index = num_batch;
			double cost = 0;
			// init
			memset2(tr.dW_1, num_e_2, num_e_1);
			memset2(tr.dW_2, num_e_3, num_e_2);
			memset(tr.dW_3, 0, sizeof(double) * num_e_3);
			memset(tr.db_1, 0, sizeof(double) * num_e_2);
			memset(tr.db_2, 0, sizeof(double) * num_e_3);
			tr.db_3 = 0;
			// end init
			for (; index < num_batch + batch_size && index < train_datas.size(); ++index) {
				data & d = train_datas[index];
				for (int i = 0; i != num_e_2; ++i) {
					double z_1 = 0;
					for (int j = 0; j != num_e_1; ++j) {
						z_1 += tr.W_1[i][j] * d.e_1[j];
					}
					E_2[i] = tanh(z_1 + tr.b_1[i]);
				}
				for (int i = 0; i != num_e_3; ++i) {
					double z_2 = 0;
					for (int j = 0; j != num_e_2; ++j) {
						z_2 += tr.W_2[i][j] * E_2[j];
					}
					E_3[i] = ReLU(z_2 + tr.b_2[i]);
				}
				y_hat = 0;
				for (int i = 0; i != num_e_3; ++i) {
					y_hat += tr.W_3[i] * E_3[i];
				}
				y_hat += tr.b_3;
				
				////////
				delta_3 = y_hat - d.res;
				for (int i = 0; i != num_e_3; ++i)
					tr.dW_3[i] += E_3[i] * delta_3;
				tr.db_3 += delta_3;
				cost += delta_3 * delta_3;

				for (int i = 0; i != num_e_3; ++i)
					delta_2[i] = delta_3 * tr.W_3[i] * dReLU(E_3[i]);
				for (int i = 0; i != num_e_3; ++i)
					for (int j = 0; j != num_e_2; ++j)
						tr.dW_2[i][j] += delta_2[i] * E_2[j];
				for (int i = 0; i != num_e_3; ++i)
					tr.db_2[i] += delta_2[i];

				for (int i = 0; i != num_e_2; ++i) {
					double d_1 = 0;
					for (int j = 0; j != num_e_3; ++j)
						d_1 += delta_2[j] * tr.W_2[j][i] * (1 - E_2[i] * E_2[i]);
					delta_1[i] = d_1;
				}

				for (int i = 0; i != num_e_2; ++i)
					for (int j = 0; j != num_e_1; ++j) {
						tr.dW_1[i][j] += delta_1[i] * d.e_1[j];
					}
				for (int i = 0; i != num_e_2; ++i)
					tr.db_1[i] += delta_1[i];
				
			}
			unsigned int N = index - num_batch;
			double J = 0;
			// J += lambda * W^2
			for (int i = 0; i != num_e_2; ++i)
				for (int j = 0; j != num_e_1; ++j)
					J += tr.W_1[i][j] * tr.W_1[i][j];
			for (int i = 0; i != num_e_3; ++i)
				for (int j = 0; j != num_e_2; ++j)
					J += tr.W_2[i][j] * tr.W_2[i][j];
			for (int i = 0; i != num_e_3; ++i)
				J += tr.W_3[i] * tr.W_3[i];
			// 0.5 * sigma [y_hat - y] ^ 2
			J = cost / (2 * N) + lambda * J / 2;
			std::cout << J << std::endl;
			for (int i = 0; i != num_e_3; ++i)
				tr.dW_3[i] = tr.dW_3[i] / N + lambda * tr.W_3[i];
			tr.db_3 /= N;

			for (int i = 0; i != num_e_3; ++i)
				for (int j = 0; j != num_e_2; ++j)
					tr.dW_2[i][j] = tr.dW_2[i][j] / N + lambda * tr.W_2[i][j];
			for (int i = 0; i != num_e_3; ++i)
				tr.db_2[i] /= N;

			for (int i = 0; i != num_e_2; ++i)
				for (int j = 0; j != num_e_1; ++j)
					tr.dW_1[i][j] = tr.dW_1[i][j] / N + lambda * tr.W_1[i][j];
			for (int i = 0; i != num_e_2; ++i)
				tr.db_1[i] /= N;
			/////
			for (int i = 0; i != num_e_2; ++i)
				for (int j = 0; j != num_e_1; ++j) {
					m_1[i][j] = beta_1 * m_1[i][j] + (1 - beta_1) * tr.dW_1[i][j];
					v_1[i][j] = beta_2 * v_1[i][j] + (1 - beta_2) * tr.dW_1[i][j] * tr.dW_1[i][j];
				}
			for (int i = 0; i != num_e_3; ++i)
				for (int j = 0; j != num_e_2; ++j) {
					m_2[i][j] = beta_1 * m_1[i][j] + (1 - beta_1) * tr.dW_2[i][j];
					v_2[i][j] = beta_2 * v_2[i][j] + (1 - beta_2) * tr.dW_2[i][j] * tr.dW_2[i][j];
				}
			for (int i = 0; i != num_e_3; ++i) {
				m_3[i] = beta_1 * m_3[i] + (1 - beta_1) * tr.dW_3[i];
				v_3[i] = beta_2 * v_3[i] + (1 - beta_2) * tr.dW_3[i] * tr.dW_3[i];
			}
			for (int i = 0; i != num_e_2; ++i) {
				m_b_1[i] = beta_1 * m_b_1[i] + (1 - beta_1) * tr.db_1[i];
				v_b_1[i] = beta_2 * v_b_1[i] + (1 - beta_2) * tr.db_1[i] * tr.db_1[i];
			}
			for (int i = 0; i != num_e_3; ++i) {
				m_b_2[i] = beta_1 * m_b_2[i] + (1 - beta_1) * tr.db_2[i];
				v_b_2[i] = beta_2 * v_b_2[i] + (1 - beta_2) * tr.db_2[i] * tr.db_2[i];
			}
			m_b_3 = beta_1 * m_b_3 + (1 - beta_1) * tr.db_3;
			v_b_3 = beta_2 * v_b_3 + (1 - beta_2) * tr.db_3 * tr.db_3;

			for (int i = 0; i != num_e_2; ++i)
				for (int j = 0; j != num_e_1; ++j) {
					tr.W_1[i][j] -= alpha * m_1[i][j] / (1 - beta_1_t) / (sqrt(v_1[i][j] / (1 - beta_2_t)) + epsilon);
				}
			for (int i = 0; i != num_e_3; ++i)
				for (int j = 0; j != num_e_2; ++j)
					tr.W_2[i][j] -= alpha * m_2[i][j] / (1 - beta_1_t) / (sqrt(v_2[i][j] / (1 - beta_2_t)) + epsilon);
			for (int i = 0; i != num_e_3; ++i)
				tr.W_3[i] -= alpha * m_3[i] / (1 - beta_1_t) / (sqrt(v_3[i] / (1 - beta_2_t)) + epsilon);
			for (int i = 0; i != num_e_2; ++i)
				tr.b_1[i] -= alpha * m_b_1[i] / (1 - beta_1_t) / (sqrt(v_b_1[i] / (1 - beta_2_t)) + epsilon);
			for (int i = 0; i != num_e_3; ++i) {
				tr.b_2[i] -= alpha * m_b_2[i] / (1 - beta_1_t) / (sqrt(v_b_2[i] / (1 - beta_2_t)) + epsilon);
			}
			tr.b_3 -= alpha * m_b_3 / (1 - beta_1_t) / (sqrt(v_b_3 / (1 - beta_2_t)) + epsilon);

			beta_1_t *= beta_1;
			beta_2_t *= beta_2;
		}
	}

	std::cout << "calcuting R_square" << std::endl;
	std::cout << "R_square of train_datas: " << getR_square(train_datas, ft.mu_res) << std::endl;
	double mean_Y = 0;
	int cnt = 1;
	for (auto &v : test_datas) {
		mean_Y += (v.res - mean_Y) / cnt;
		++cnt;
	}
	std::cout << "R_square of test_datas: " << getR_square(test_datas, mean_Y) << std::endl;

	delete[] m_1;
	delete[] m_2;
	delete[] m_3;
	delete[] m_b_1;
	delete[] m_b_2;
	delete[] v_1;
	delete[] v_2;
	delete[] v_3;
	delete[] v_b_1;
	delete[] v_b_2;
	delete[] E_2;
	delete[] E_3;
	delete[] delta_1;
	delete[] delta_2;

}

double DataTrainer::getR_square(std::vector<data> & vec, double mean_Y) {
	double *E_2 = new double[num_e_2];
	double *E_3 = new double[num_e_3];
	//E_2[NUM_E_2], E_3[NUM_E_3];
	double y_hat;
	double dY_hat = 0, dY = 0;
	for (auto &d : vec) {
		for (int i = 0; i != num_e_2; ++i) {
			double z_1 = 0;
			for (int j = 0; j != num_e_1; ++j) {
				z_1 += tr.W_1[i][j] * d.e_1[j];
			}
			E_2[i] = tanh(z_1 + tr.b_1[i]);
		}
		for (int i = 0; i != num_e_3; ++i) {
			double z_2 = 0;
			for (int j = 0; j != num_e_2; ++j) {
				z_2 += tr.W_2[i][j] * E_2[j];
			}
			E_3[i] = ReLU(z_2 + tr.b_2[i]);
		}
		y_hat = 0;
		for (int i = 0; i != num_e_3; ++i) {
			y_hat += tr.W_3[i] * E_3[i];
		}
		y_hat += tr.b_3;
		dY_hat += square(d.res - y_hat);
		dY += square(d.res - mean_Y);
		//d.res -= mean_Y;
	}
	delete[] E_2;
	delete[] E_3;
	return 1 - dY_hat / dY;
}

void DataTrainer::transform(std::vector<data> & vec, interval & itv) {
	itv.tmax = new double[pre_num];
	itv.tmin = new double[pre_num];
	itv.pov = new double[pre_num_res * (pre_num_res - 1) / 2];
	for (int i = 0; i != pre_num; ++i) {
		itv.tmax[i] = -11111111.0;
		itv.tmin[i] = 11111111.0;
	}
	for (auto & v : vec) {
		for (size_t i = 0; i != num_e_1; ++i) {
			v.e_1[i] = (v.e_1[i] - ft.mu[i]) / ft.sigma[i];
			if (i >= pre_num_res) continue;
			if (v.e_1[i] > itv.tmax[i]) itv.tmax[i] = v.e_1[i];
			if (v.e_1[i] < itv.tmin[i]) itv.tmin[i] = v.e_1[i];
		}
		//v.res = (v.res - ft.mu_res) / ft.sigma_res;
		if (v.res > itv.tmax[pre_num_res]) itv.tmax[pre_num_res] = v.res;
		if (v.res < itv.tmin[pre_num_res]) itv.tmin[pre_num_res] = v.res;
	}
	/*
	for (auto &v : vec) {
		v.res += 1.5 - itv.tmin[PRE_NUM_RES];
	}
	*/
}

void DataTrainer::show() {
	
	//
	AW = new TR::Axis::AxisWindow(10, 10, 800, 620, "Data Graph", train_datas, train_itv, name_list, res_index);
	AW->init();
	Fl::run();
}

void TR::read(FileReader * fr, DataTrainer * dt) {
	std::cout << "Loading datas ..." << std::endl;
	fr->read(dt->train_datas, dt->test_datas, dt->ft, dt->res_index, dt->name_list);
	// transform
	dt->transform(dt->train_datas, dt->train_itv);
	dt->transform(dt->test_datas, dt->test_itv);
	interval &itv = dt->train_itv;
	for (int i = 0; i != pre_num_res * (pre_num_res - 1) / 2; ++i) {
		itv.pov[i] = 0.0;
	}
	unsigned int cnt = 1;
	for (auto &v : dt->train_datas) {
		for (size_t i = 0; i != pre_num_res; ++i)
			for (size_t j = i + 1; j < pre_num_res; ++j) {
				int index = i * (pre_num_res - 1) - i * (i + 1) / 2 + j - 1;
				itv.pov[index] += (v.e_1[i] * v.e_1[j] - itv.pov[index]) / cnt;
			}
		++cnt;
	}
	std::cout << "Loading finished." << std::endl;
}
