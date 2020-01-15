#define MM
#ifdef MM
#include "Trainer.h"
#include "Drawer.h"
#include <iostream>
using namespace std;
int main() {

	// 用于预处理数据
	cout << "input filename:";
	string filename;// = "Asteroid_Updated.csv";
	getline(cin, filename);
	cout << "pre-creating datas...";
	TR::FileReader *fr = new TR::FileReader(filename, "diameter", 30000, 1);
	/*
		for test
	*/
	//std::string filename = "test.csv";
	//std::string filename = "Asteroid_train.csv";
	//TR::FileReader *fr = new TR::FileReader(filename);
	TR::DataTrainer *dt = new TR::DataTrainer();
	TR::read(fr, dt);
	// datas : use 38 MB memories
	// epoch = 100, batch_size = 256, alpha = 0.005, lambda = 0
	dt->traindatas(100, 256, 0.005, 0);

	dt->show();
	

	return 0;
}

#endif
