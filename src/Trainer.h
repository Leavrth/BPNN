#pragma once

#include "DataStruct.h"
#include "Drawer.h"
namespace TR {
	class FileReader {
	public:
		explicit FileReader(const std::string & str, const std::string & res_str, size_t line_num, size_t np = 10);
		void read(std::vector<data> &, std::vector<data> & datas2, fitter & ft, size_t &res_index, std::map<std::string, size_t> &name_list);
		~FileReader();
	private:
		int hashcode(const std::string &str);
		void del(size_t res_index, size_t iter_cur_num);
		size_t num_of_partitions;
		std::string FILE_RD;
		std::set<std::string> FILE_NAMES;
		std::map<std::string, size_t> name_list; // 名字与索引的map,用于下面数据结构
		size_t *isnan_cnt;
		std::string res_str;
	};
	class DataTrainer {
	public:
		friend void read(FileReader *, DataTrainer *);
		void traindatas(int epoch, int batch_size, double alpha, double lambda);
		void show();
	private:
		void transform(std::vector<data> & vec, interval & itv);
		double getR_square(std::vector<data> & vec, double mean_Y);
		std::vector<data> train_datas;
		std::vector<data> test_datas;
		fitter ft;
		interval train_itv, test_itv;
		trainres tr;
		size_t res_index;
		std::map<std::string, size_t> name_list;
		
		TR::Axis::AxisWindow * AW;
	};
	
	void read(FileReader * fr, DataTrainer * dt);
	
}