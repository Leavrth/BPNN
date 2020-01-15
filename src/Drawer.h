#pragma once
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include "DataStruct.h"

namespace TR {
	
	namespace Axis {
		class  AxisWindow : public Fl_Window {
		public:
			AxisWindow(int x, int y, int w, int h, const std::string & titile, const std::vector<TR::data> &, const TR::interval &, const std::map<std::string, size_t> &, size_t);
			void init();
			virtual ~AxisWindow() { delete[] density; delete[] stitle; }
			void nextG() { 
				graph_cur = (graph_cur + 1) % 3; 
				if (graph_cur == 1) nextPt->set_visible(); 
				else nextPt->clear_visible();
				if (graph_cur == 0) this->resize(xxx, yyy, width, height + 200);
				else this->resize(xxx, yyy, width, height);
			}
			void nextP() { parameter_cur = (parameter_cur + 1) % pre_num_res; }
		protected:
			void draw();
		private:
			void drawAxis(int x, int y, int w, int h, const TR::Axisdata & Axsd);
			std::vector<TR::Axisdata> Axsds;
			const std::vector<TR::data> & ref_datas;
			const TR::interval & ref_itv;
			const std::map<std::string, size_t> & name_list;
			size_t parameter_cur;
			size_t res_index;
			int xxx, yyy;
			int width, height;
			double **density; //density[PRE_NUM][AXIS_DIV];
			Fl_Button *nextBt;
			Fl_Button *nextPt;
			int graph_cur;
			std::string *parameter_name;
			std::string *stitle; //stitle[PRE_NUM];
		};
	}
}