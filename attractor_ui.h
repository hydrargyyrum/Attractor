#pragma once
#include <vector>
using namespace std;

struct attractor_show_args {
	bool example_to_show = true;
	bool example_continiously_draw = true;
	vector<vector<double>> points;
};

int attractor_show(attractor_show_args args);