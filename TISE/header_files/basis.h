#pragma once

#include <iostream>
#include <vector>

namespace basis 
{
	double B(int i, int order, const std::vector<double>& knots, double x);
	double dB(int i, int order, const std::vector<double>& knots, double x);
	std::vector<double> linear_knots(int N, int k, double xmax);
	double overlap_matrix_element(int i, int j, int degree, const std::vector<double>& knots);
}
