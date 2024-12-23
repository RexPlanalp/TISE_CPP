#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <functional>

namespace basis 
{	
	std::vector<std::complex<double>> linear_knots(int n_basis, int degree, double rmax);
	std::complex<double> B(int i, int degree, const std::vector<std::complex<double>>& knots, std::complex<double> x);
	std::complex<double> dB(int i, int degree, const std::vector<std::complex<double>>& knots, std::complex<double> x);
	std::complex<double> compute_matrix_element(int i, int j, int degree, const std::vector<std::complex<double>>& knots,std::function<std::complex<double>(int, int, const std::vector<std::complex<double>>&, std::complex<double>)> integrand);
	std::complex<double> overlap_matrix_element(int i, int j, int degree, const std::vector<std::complex<double>>& knots);
	std::complex<double> kinetic_matrix_element(int i, int j, int degree, const std::vector<std::complex<double>>& knots);
	std::complex<double> inverse_r2_matrix_element(int i, int j, int degree, const std::vector<std::complex<double>>& knots);
	std::complex<double> inverse_r_matrix_element(int i, int j, int degree, const std::vector<std::complex<double>>& knots);
	int save_bsplinee_basis(int n_basis, int degree, const std::vector<std::complex<double>>& knots, int Nx, double xmax);
}
