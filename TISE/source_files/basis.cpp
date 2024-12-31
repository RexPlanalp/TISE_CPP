#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <functional>
#include <complex>
#include <unordered_map>
#include <petscsys.h>



namespace basis
{

	std::array<double,2> roots_two = {-0.57735027,0.57735027};
	std::array<double,2> weights_two = {1, 1};

	std::array<double,3> roots_three = {-0.77459667, 0.0, 0.77459667};
	std::array<double,3> weights_three = {0.55555556, 0.88888889, 0.55555556};

	std::array<double,4> roots_four = {-0.86113631, -0.33998104,  0.33998104,  0.86113631};
	std::array<double,4> weights_four = {0.34785485, 0.65214515, 0.65214515, 0.34785485};

	std::array<double,5> roots_five = {-0.90617985, -0.53846931,  0.0, 0.53846931, 0.90617985};
	std::array<double,5> weights_five = {0.23692689, 0.47862867, 0.56888889, 0.47862867, 0.23692689};

	std::array<double,6> roots_six = {-0.93246951, -0.66120939, -0.23861919,  0.23861919,  0.66120939,  0.93246951};
	std::array<double,6> weights_six = {0.17132449, 0.36076157, 0.46791393, 0.46791393, 0.36076157, 0.17132449};

	std::array<double, 7> roots_seven = {-0.94910791,-0.74153119,-0.40584515,0,0.40584515,0.74153119,0.94910791};
	std::array<double, 7> weights_seven = { 0.12948497,0.27970539,0.38183005,0.41795918,0.38183005,0.27970539,0.12948497};

	std::array<double, 8> roots_eight = { -0.96028986, -0.79666648, -0.52553241, -0.18343464, 0.18343464, 0.52553241, 0.79666648, 0.96028986};
	std::array<double, 8> weights_eight = {0.10122854, 0.22238103, 0.31370665, 0.36268378, 0.36268378, 0.31370665, 0.22238103, 0.10122854};


	std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>> gauss_quadrature =
	{
		{1, {std::vector<double>{basis::roots_two.begin(), basis::roots_two.end()}, 
			std::vector<double>{basis::weights_two.begin(), basis::weights_two.end()}}},
		{2, {std::vector<double>{basis::roots_three.begin(), basis::roots_three.end()}, 
			std::vector<double>{basis::weights_three.begin(), basis::weights_three.end()}}},
		{3, {std::vector<double>{basis::roots_four.begin(), basis::roots_four.end()}, 
			std::vector<double>{basis::weights_four.begin(), basis::weights_four.end()}}},
		{4, {std::vector<double>{basis::roots_five.begin(), basis::roots_five.end()}, 
			std::vector<double>{basis::weights_five.begin(), basis::weights_five.end()}}},
		{5, {std::vector<double>{basis::roots_six.begin(), basis::roots_six.end()}, 
			std::vector<double>{basis::weights_six.begin(), basis::weights_six.end()}}},
		{6, {std::vector<double>{basis::roots_seven.begin(), basis::roots_seven.end()}, 
			std::vector<double>{basis::weights_seven.begin(), basis::weights_seven.end()}}},
		{7, {std::vector<double>{basis::roots_eight.begin(), basis::roots_eight.end()}, 
			std::vector<double>{basis::weights_eight.begin(), basis::weights_eight.end()}}}
	};

	PetscScalar R_x(PetscScalar x, double R0, double eta)
	{
		if (eta == 0)
		{
			return x;
		}

		if (x.real() < R0)
		{
			return x;
		}
		else
		{
			return R0 + (x - R0) * std::exp(PetscScalar(0, M_PI * eta));
		}
	}

	std::vector<PetscScalar> R_k(const std::vector<PetscScalar>& knots, double R0, double eta)
	{
		if (eta == 0)
		{
			return knots;
		}

		std::vector<PetscScalar> modified(knots.size());
		std::transform(knots.begin(), knots.end(), modified.begin(),
					[R0, eta](const PetscScalar& x) {
						return R_x(x, R0, eta);
					});
		return modified;
	}

	PetscScalar R_w(PetscScalar x, double w, double R0, double eta)
	{

		if (eta == 0)
		{
			return w;
		}


		if (x.real() < R0)
		{
			return w;
		}
		else
		{
			return w * std::exp(PetscScalar(0, M_PI * eta));			
		}
	}

	std::vector<PetscScalar> linear_knots(int n_basis, int degree, double rmax) 
	{	
		int order = degree + 1;


		std::vector<PetscScalar> knots;
		int N_knots = n_basis + order;

		int N_middle = N_knots - 2 * (order - 2);
		double step_size = rmax / (N_middle-1);
		std::vector<PetscScalar> knots_middle;
		for (int idx = 0; idx < N_middle; ++idx) 
		{
			knots_middle.push_back(idx * step_size);
		}
		knots_middle.back() = rmax;



		std::vector<PetscScalar> knots_start(order - 2, 0.0);
		std::vector<PetscScalar> knots_end(order - 2, rmax);

		knots.insert(knots.end(), knots_start.begin(), knots_start.end());
		knots.insert(knots.end(), knots_middle.begin(), knots_middle.end());
		knots.insert(knots.end(), knots_end.begin(), knots_end.end());
		return knots;
	}

	PetscScalar B(PetscInt i, PetscInt degree,  const std::vector<PetscScalar>& knots, PetscScalar x)
	{
		if (degree == 0)
		{
			return (knots[i].real() <= x.real() && x.real() < knots[i + 1].real()) ? 1.0 : 0.0;
		}

		PetscScalar denom1 = knots[i + degree] - knots[i];
		PetscScalar denom2 = knots[i + degree + 1] - knots[i + 1];

		PetscScalar term1 = 0.0;
		PetscScalar term2 = 0.0;

		if (denom1.real() > 0)
		{
			term1 = (x - knots[i]) / (denom1) * B(i, degree - 1, knots, x);
		}

		if (denom2.real() > 0)
		{
			term2 = (knots[i + degree + 1] - x) / (denom2)  * B(i + 1, degree - 1, knots, x);
		}

		return term1+term2;
	}

	PetscScalar dB(int i, int degree, const std::vector<PetscScalar>& knots, PetscScalar x)
	{
		if (degree == 0)
		{
			return 0.0;
		}

		PetscScalar denom1 = knots[i + degree] - knots[i];
		PetscScalar denom2 = knots[i + degree + 1] - knots[i + 1];

		PetscScalar term1 = 0.0;
		PetscScalar term2 = 0.0;

		if (denom1.real() > 0)
		{
			term1 = (PetscScalar(degree)) / (denom1) * B(i, degree - 1, knots, x);
		}

		if (denom2.real() > 0)
		{
			term2 = (-PetscScalar(degree)) / (denom2)*B(i + 1, degree - 1, knots, x);
		}

		return term1 + term2;
	}

	PetscScalar compute_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots, std::function<PetscScalar(int, int, const std::vector<PetscScalar>&, PetscScalar)> integrand, double R0, double eta) // Note eta here
	{
		
		std::unordered_map<int, std::pair<std::vector<double>, std::vector<double>>>::iterator it = gauss_quadrature.find(degree);
		std::vector<double>& roots = it->second.first;
		std::vector<double>& weights = it->second.second;

		std::vector<PetscScalar> knots_modified = R_k(knots, R0, eta);

		PetscScalar total = 0.0;

		int lower = std::min(i, j);
		int upper = std::max(i, j);

		for (int k = lower; k <= upper + degree; ++k) {
			PetscScalar a = knots[k];    
			PetscScalar b = knots[k + 1]; 

			if (a == b)
				continue;

			for (size_t r = 0; r < roots.size(); ++r) {
				PetscScalar xi = 0.5 * (b - a) * roots[r] + 0.5 * (b + a); 
				double weight = weights[r];

				
				PetscScalar xi_modified = R_x(xi, R0, eta);
				PetscScalar weight_modified = R_w(xi,weight,R0,eta);
				total += weight_modified * integrand(i, j, knots_modified, xi_modified) * (b - a) * 0.5;
			}
		}


		return total;
	}

	PetscScalar overlap_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots, double R0, double eta) 
	{
		return compute_matrix_element(
			i, 
			j, 
			degree, 
			knots, 
			[degree, R0, eta](int i, int j, const std::vector<PetscScalar>& knots, PetscScalar xi)
			{
			PetscScalar Bi = basis::B(i, degree, knots, xi); PetscScalar Bj = basis::B(j, degree, knots, xi);
			return Bi * Bj;
			}, 
			R0, 
			eta // Pass R0 and w explicitly here
		);
	}

	PetscScalar kinetic_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots,double R0,double eta) {
		return compute_matrix_element(i, j, degree, knots, 
			[degree](int i, int j, const std::vector<PetscScalar>& knots, PetscScalar xi) {
				PetscScalar Bi = basis::dB(i, degree, knots, xi);
				PetscScalar Bj = basis::dB(j, degree, knots, xi);
				return 0.5 * Bi * Bj;
			},R0,eta);
	}

	PetscScalar inverse_r2_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots,double R0,double eta) {
		return compute_matrix_element(i, j, degree, knots, 
			[degree](int i, int j, const std::vector<PetscScalar>& knots, PetscScalar xi) {
				PetscScalar Bi = basis::B(i, degree, knots, xi);
				PetscScalar Bj = basis::B(j, degree, knots, xi);
				return Bi * Bj / (xi * xi + 1E-25);
			},R0,eta);
	}

	PetscScalar inverse_r_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots,double R0,double eta) {
		return compute_matrix_element(i, j, degree, knots, 
			[degree](int i, int j, const std::vector<PetscScalar>& knots, PetscScalar xi) {
				PetscScalar Bi = basis::B(i, degree, knots, xi);
				PetscScalar Bj = basis::B(j, degree, knots, xi);
				return Bi * Bj / (xi + 1E-25);
			},R0,eta);
	}

	int save_bsplinee_basis(int n_basis, int degree, const std::vector<PetscScalar>& knots, int Nx, double xmax,double R0,double eta)
	{
		std::vector<PetscScalar> x_vector;
		double step_size = xmax / (Nx - 1);

		// Generate x_vector
		for (int idx = 0; idx < Nx; ++idx)
		{
			x_vector.push_back(idx * step_size);
		}

		// Open the file for writing
		std::ofstream file("bsplines.txt");
		if (!file.is_open())
		{
			std::cerr << "Error: Could not open file for writing." << std::endl;
			return -1;
		}

		// Loop through each basis function
		for (int i = 0; i < n_basis; ++i)
		{
			for (int idx = 0; idx < Nx; ++idx)
			{

				std::vector<PetscScalar> knots_modified = R_k(knots,R0,eta);
				PetscScalar xi_modified = R_x(x_vector[idx],R0,eta);
				// Compute the B-spline value for the given x
				PetscScalar y = B(i, degree, knots_modified, xi_modified);

				// Write the real and imaginary parts as columns
				file << y.real() << "\t" << y.imag() << "\n";
			}
			// Add a blank line between basis functions for better readability (optional)
			file << "\n";
		}

		// Close the file
		file.close();
		return 0;
	}
}