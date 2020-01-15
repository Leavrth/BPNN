#include "Drawer.h"
using namespace TR::Axis;

AxisWindow::AxisWindow(int x, int y, int w, int h, const std::string & title, const std::vector<TR::data> & datas, const TR::interval & itv, const std::map<std::string, size_t> & name_list, size_t res_index)
	: Fl_Window(x, y, w, h, title.c_str()), xxx(x), yyy(y), width(w), height(h), graph_cur(0), ref_datas(datas), ref_itv(itv), name_list(name_list), res_index(res_index), parameter_cur(0) {
	density = new double*[pre_num];
	for (int i = 0; i != pre_num; ++i) density[i] = new double[AXIS_DIV];
	stitle = new std::string[pre_num];
	parameter_name = new std::string[pre_num_res];
	begin();
	nextBt = new Fl_Button(w-20, h - 20, 20, 20, "@->");
	nextPt = new Fl_Button(w - 20, 0, 20, 20, "@->");
	nextPt->clear_visible();
	nextBt->callback([](Fl_Widget *o, void *v) {
		AxisWindow * AW = static_cast<AxisWindow *>(v);
		AW->nextG();
		AW->redraw();
	}, this);
	nextPt->callback([](Fl_Widget *o, void *v) {
		AxisWindow * AW = static_cast<AxisWindow *>(v);
		AW->nextP();
		AW->redraw();
	}, this);
	
	this->color(fl_rgb_color(96, 143, 159));
	nextBt->color(fl_rgb_color(39, 72, 98));
	nextBt->color2(fl_rgb_color(39, 72, 98));
	nextBt->down_color(fl_rgb_color(39, 72, 98));
	nextBt->labelcolor(fl_rgb_color(34, 54, 32));
	nextPt->color(fl_rgb_color(39, 72, 98));
	nextPt->color2(fl_rgb_color(39, 72, 98));
	nextPt->down_color(fl_rgb_color(39, 72, 98));
	nextPt->labelcolor(fl_rgb_color(34, 54, 32));
	//draw();
	end();
	resizable(this);
	show();
	
}


void AxisWindow::init() {
	size_t **cnt = new size_t*[pre_num];
	for (int i = 0; i != pre_num; ++i) {
		cnt[i] = new size_t[AXIS_DIV]();
	}
	for (auto & p : name_list) {
		stitle[p.second] = std::string("Density of ").append(p.first);
		if (p.second != res_index) {
			size_t cor = p.second > res_index ? p.second - 1 : p.second;
			parameter_name[cor] = p.first;
		}	
	}

	for (const auto &v : ref_datas) {
		for (int i = 0; i != pre_num_res; ++i) {
			int j = static_cast<int>((v.e_1[i] - ref_itv.tmin[i]) / (ref_itv.tmax[i] - ref_itv.tmin[i]) * AXIS_DIV);
			j = j == AXIS_DIV ? j - 1 : j;
			++cnt[i][j];
		}
		int j = static_cast<int>((v.res - ref_itv.tmin[pre_num_res]) / (ref_itv.tmax[pre_num_res] - ref_itv.tmin[pre_num_res]) * AXIS_DIV);
		j = j == AXIS_DIV ? j - 1 : j;
		++cnt[pre_num_res][j];
	}

	double *y_min = new double[pre_num];
	double *y_max = new double[pre_num];
	for (int i = 0; i != pre_num; ++i) {
		y_min[i] = 11111111.0;
		y_max[i] = -11111111.0;
	}
	for (int i = 0; i != pre_num; ++i)
		for (int j = 0; j != AXIS_DIV; ++j) {
			density[i][j] = 1.0 * cnt[i][j] / ref_datas.size();
			if (density[i][j] > y_max[i]) y_max[i] = density[i][j];
			if (density[i][j] < y_min[i]) y_min[i] = density[i][j];
		}

	for (int i = 0; i != pre_num; ++i)
		Axsds.emplace_back(Axisdata(stitle[i].c_str(), "x", "density",
			ref_itv.tmin[i], ref_itv.tmax[i], y_min[i], y_max[i],
			AXIS_DIV, AXIS_DIV, density[i]));
	


	this->redraw();
	this->resize(xxx, yyy, width, height + 200);
}

void AxisWindow::draw() {
	Fl_Window::draw();
	Fl_Color oldc = fl_color();
	// there is no good portable way of retrieving the current style
	fl_color(fl_rgb_color(254, 248, 230));            // set color
	fl_line_style(FL_SOLID, 1); // set style
	/*
	for (int i = 0; i != width; ++i)
		fl_line(i, 0, i, height-20);
		*/
	fl_color(FL_BLACK);
	switch (graph_cur) {
	case 0: {
		this->label("Density of parameters");
		int c = 0, l = 0;  // l : 0 - 3
		for (auto & v : Axsds) {
			drawAxis(l, c, width / 4, width / 4, v);
			if (++l == 4) { l = 0; ++c; }
			}
		break;
		}
	case 1: {
		this->label("Relationship between the result and one parameter");
		fl_font(FL_COURIER_BOLD, 14);
		fl_draw("Relationship between the result and one parameter", 200, 24);
		fl_font(FL_COURIER, 12);
		fl_draw("diameter", 10, 10);
		fl_draw(parameter_name[parameter_cur].c_str(), width / 2 - 50, height - 10);
		fl_line(20, 20, 20, height - 15);
		fl_line(15, height - 20, width - 20, height - 20);
		fl_line(20, 20, 15, 29);
		fl_line(20, 20, 25, 29);
		fl_line(width - 20, height - 20, width - 29, height - 15);
		fl_line(width - 20, height - 20, width - 29, height - 25);
		double ymin = ref_itv.tmin[pre_num_res], ydelta = ref_itv.tmax[pre_num_res] - ymin;
		double xmin = ref_itv.tmin[parameter_cur], xdelta = ref_itv.tmax[parameter_cur] - xmin;
		const int Y_itv = height - 45;
		const int X_itv = width - 40;
		for (auto &v : ref_datas) {
			double x = 25 + (v.e_1[parameter_cur] - xmin) / xdelta * X_itv;
			double y = height - 20 - (v.res - ymin) / ydelta * Y_itv;
			fl_circle(x, y, 2);
		}
		break;
		}
	case 2: {
		this->label("Heat map of parameters(Pearson correlation coefficient)");
		fl_font(FL_COURIER, 12);
		const size_t l_ = 35;
		/*
		for (int i = 25, step = l_, endlth = i + step * pre_num_res, j = 0; j != pre_num; ++j, i += step) {
			fl_line(25, i, endlth, i);
			fl_line(i, 25, i, endlth);
		}
		*/
		static std::string *num = build_num();
		
		for (int i = 0, j = 47, step = l_, jj = 20; i != pre_num_res; ++i, j += step, jj += 14) {
			fl_draw(num[i].c_str(), 10, j + 6);
			fl_draw(num[i].c_str(), j, 20);
			fl_draw(num[i].c_str(), 580, jj);
			fl_draw(parameter_name[i].c_str(), 600, jj);
		}
		for (int i = 0; i != 400; i++) {
			uchar scolor = 255 - i * 256 / 400;
			fl_color(fl_rgb_color(scolor, scolor, scolor));
			fl_line(730, i + 30, 745, i + 30);
		}
		fl_draw("1", 750, 33);
		fl_draw("0", 750, 232);
		fl_draw("-1", 750, 430);
		for (int i = 0; i != pre_num_res; ++i) {
			fl_color(fl_rgb_color(255, 255, 255));
			for (int xx = 25 + l_ * i, k = xx + l_; k > xx; --k) {
				fl_line(xx, k, xx + l_, k);
			}
			for (int j = i + 1; j < pre_num_res; ++j) {
				int xx = 25 + l_ * i;
				int yy = 25 + l_ * j;
				int index = i * (pre_num_res - 1) - i * (i + 1) / 2 + j - 1;
				uchar scolor = (1 + ref_itv.pov[index]) * 255 / 2;
				fl_color(fl_rgb_color(scolor, scolor, scolor));
				for (int k = xx + l_; k > xx; --k) {
					fl_line(yy, k, yy + l_, k);
					fl_line(k, yy, k, yy + l_);
				}
			}
		}

		break;
		}
	}
	fl_color(oldc);      // reset color (to previous)
	fl_line_style(0);    // reset line style to default

}

void AxisWindow::drawAxis(int x, int y, int w, int h, const TR::Axisdata & Axsd) {
	int xx = x * w; int yy = y * h;
	/*
		(xx, yy)    w
			¡¤---------------¡¤
			|  ¡ü  title		|
			|  |			|
		  h |y |			|
			|  +----------¡ú	|
			|		x		|
			¡¤---------------¡¤
	*/
	fl_font(FL_COURIER_BOLD, 14);
	fl_draw(Axsd.title, xx + 30, yy + (h > 14 ? 14 : h));
	fl_font(FL_COURIER, 12);
	fl_draw("y", xx + 10, yy + h / 2);
	fl_draw("x", xx + w / 2, yy + h - 10);
	fl_line(xx + 20, yy + 20, xx + 20, yy + h - 15);
	fl_line(xx + 15, yy + h - 20, xx + w - 20, yy + h - 20);
	fl_line(xx + 20, yy + 20, xx + 15, yy + 29);
	fl_line(xx + 20, yy + 20, xx + 25, yy + 29);
	fl_line(xx + w - 20, yy + h - 20, xx + w - 29, yy + h - 15);
	fl_line(xx + w - 20, yy + h - 20, xx + w - 29, yy + h - 25);
	int x_begin = xx + 20, x_end = xx + w - 30, x_step = (w - 50) / Axsd.X_div_num;
	for (int j = 0; x_begin < x_end; x_begin += x_step,++j)
		for (int i = x_begin + x_step - 1; i >= x_begin; --i)
			fl_line(i, yy + h - 20, i, yy + h - 20 - (Axsd.Y[j] - Axsd.y_min) / (Axsd.y_max - Axsd.y_min)*(h - 50));
}