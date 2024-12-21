#include <iostream>
#include <vector>
#include "basis.h"
#include <fstream>
#include <iomanip>
#include <chrono>

int main()
{
    double xmax = 10;
    int n_basis = 30;
    int order = 7;
	int degree = order - 1;


    std::vector knots = basis::linear_knots(n_basis, order, xmax);


	int Nx = 1000;
	std::vector<double> x_vector;
	double step_size = xmax / (Nx - 1);


	for (int idx = 0; idx < Nx; ++idx)
	{
		x_vector.push_back(idx * step_size);
	}

	std::vector<double> y_vector;


	
	for (int i = 0; i < n_basis; ++i)
	{
		for (int idx = 0; idx < Nx; ++idx)
		{
			double y = basis::B(i, degree, knots, x_vector[idx]);
			y_vector.push_back(y);
		}
	}


	
	
	

	std::ofstream file("output.txt");

	for (const auto& y : y_vector)
	{
		file << y << std::endl;
	}


	auto start = std::chrono::high_resolution_clock::now();
	double integral_value = basis::overlap_matrix_element(15, 15, degree, knots);
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::fixed << std::setprecision(20) << integral_value << std::endl;
	std::chrono::duration<double, std::milli> duration = end - start;

	// Print the duration
	std::cout << "Time taken: " << duration.count() << " ms" << std::endl;
	return 0;




}

