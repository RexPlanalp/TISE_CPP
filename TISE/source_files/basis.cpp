#include <iostream>
#include <vector>
#include <array>
#include <immintrin.h>
#include <unordered_map>
#include <algorithm>

namespace basis
{

	std::array<double, 7> roots_seven = {-0.94910791,-0.74153119,-0.40584515,0,0.40584515,0.74153119,0.94910791};
	std::array<double, 7> weights_seven = { 0.12948497,0.27970539,0.38183005,0.41795918,0.38183005,0.27970539,0.12948497};

	std::array<double, 8> roots_eight = { -0.96028986, -0.79666648, -0.52553241, -0.18343464, 0.18343464, 0.52553241, 0.79666648, 0.96028986};
	std::array<double, 8> weights_eight = {0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854};

	double B(int i, int degree, const std::vector<double>& knots, double x)
	{
		if (degree == 0)
		{
			return (knots[i] <= x && x < knots[i + 1]) ? 1.0 : 0.0;
		}

		double denom1 = knots[i + degree] - knots[i];
		double denom2 = knots[i + degree + 1] - knots[i + 1];

		double term1 = 0.0;
		double term2 = 0.0;

		if (denom1 != 0)
		{
			term1 = (x - knots[i]) / (denom1) * B(i, degree - 1, knots, x);
		}

		if (denom2 != 0)
		{
			term2 = (knots[i + degree + 1] - x) / (denom2)  *B(i + 1, degree - 1, knots, x);
		}

		return term1+term2;
	}

	double dB(int i, int degree, const std::vector<double>& knots, double x)
	{
		if (degree == 0)
		{
			return 0.0;
		}

		double denom1 = knots[i + degree] - knots[i];
		double denom2 = knots[i + degree + 1] - knots[i + 1];

		double term1 = 0.0;
		double term2 = 0.0;

		if (denom1 != 0)
		{
			term1 = (degree) / (denom1) * B(i, degree - 1, knots, x);
		}

		if (denom2 != 0)
		{
			term2 = (-degree) / (denom2)*B(i + 1, degree - 1, knots, x);
		}

		return term1 + term2;
	}

	std::vector<double> linear_knots(int n_basis, int degree, double rmax) 
	{	
		int order = degree + 1;


		std::vector<double> knots;
		int N_knots = n_basis + order;

		int N_middle = N_knots - 2 * (order - 2);
		double step_size = rmax / (N_middle-1);
		std::vector<double> knots_middle;
		for (int idx = 0; idx < N_middle; ++idx) 
		{
			knots_middle.push_back(idx * step_size);
		}
		knots_middle.back() = rmax;



		std::vector<double> knots_start(order - 2, 0.0);
		std::vector<double> knots_end(order - 2, rmax);

		knots.insert(knots.end(), knots_start.begin(), knots_start.end());
		knots.insert(knots.end(), knots_middle.begin(), knots_middle.end());
		knots.insert(knots.end(), knots_end.begin(), knots_end.end());



		return knots;
	}

	double overlap_matrix_element(int i, int j, int degree, const std::vector<double>& knots)
	{
		double total = 0.0;

		int lower = std::min(i, j);
		int upper = std::max(i, j);

		for (int k = lower; k <= upper + degree; ++k)
		{
			double a = knots[k];    
			double b = knots[k + 1]; 

			if (a == b)
				continue;

			for (size_t r = 0; r < basis::roots_seven.size(); ++r)
			{
				double xi = 0.5 * (b - a) * basis::roots_seven[r] + 0.5 * (b + a); 
				double weight = basis::weights_seven[r];

				double Bi = basis::B(i, degree, knots, xi);
				double Bj = basis::B(j, degree, knots, xi);

				total += weight * Bi * Bj * (b - a) * 0.5;
				
			}
		}

		return total;
	}

	double kinetic_matrix_element(int i, int j, int degree, const std::vector<double>& knots)
	{
		double total = 0.0;

		int lower = std::min(i, j);
		int upper = std::max(i, j);

		for (int k = lower; k <= upper + degree; ++k)
		{
			double a = knots[k];    
			double b = knots[k + 1]; 

			if (a == b)
				continue;

			for (size_t r = 0; r < basis::roots_seven.size(); ++r)
			{
				double xi = 0.5 * (b - a) * basis::roots_seven[r] + 0.5 * (b + a); 
				double weight = basis::weights_seven[r];

				double Bi = basis::dB(i, degree, knots, xi);
				double Bj = basis::dB(j, degree, knots, xi);

				total += weight * Bi * Bj * (b - a) * 0.5;
				
			}
		}

		return total * 0.5;
	}

	double inverse_r2_matrix_element(int i, int j, int degree, const std::vector<double>& knots)
	{
		double total = 0.0;

		int lower = std::min(i, j);
		int upper = std::max(i, j);

		for (int k = lower; k <= upper + degree; ++k)
		{
			double a = knots[k];    
			double b = knots[k + 1]; 

			if (a == b)
				continue;

			for (size_t r = 0; r < basis::roots_seven.size(); ++r)
			{
				double xi = 0.5 * (b - a) * basis::roots_seven[r] + 0.5 * (b + a); 
				double weight = basis::weights_seven[r];

				double Bi = basis::B(i, degree, knots, xi);
				double Bj = basis::B(j, degree, knots, xi);

				total += weight * Bi * Bj * (b - a) * 0.5 / (xi*xi + 1E-25);
				
			}
		}

		return total;
	}

	double inverse_r_matrix_element(int i, int j, int degree, const std::vector<double>& knots)
	{
		double total = 0.0;

		int lower = std::min(i, j);
		int upper = std::max(i, j);

		for (int k = lower; k <= upper + degree; ++k)
		{
			double a = knots[k];    
			double b = knots[k + 1]; 

			if (a == b)
				continue;

			for (size_t r = 0; r < basis::roots_seven.size(); ++r)
			{
				double xi = 0.5 * (b - a) * basis::roots_seven[r] + 0.5 * (b + a); 
				double weight = basis::weights_seven[r];

				double Bi = basis::B(i, degree, knots, xi);
				double Bj = basis::B(j, degree, knots, xi);

				total += weight * Bi * Bj * (b - a) * 0.5 / (xi + 1E-25);
				
			}
		}

		return total;
	}

}