#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <functional>

namespace basis 
{	
	std::vector<PetscScalar> linear_knots(int n_basis, int degree, double rmax);
	PetscScalar B(int i, int degree, const std::vector<PetscScalar>& knots, PetscScalar x);
	PetscScalar dB(int i, int degree, const std::vector<PetscScalar>& knots, PetscScalar x);


	PetscScalar compute_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots, std::function<PetscScalar(int, int, const std::vector<PetscScalar>&, PetscScalar)> integrand, double R0, double eta);
	PetscScalar overlap_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots, double R0, double w);
	PetscScalar kinetic_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots,double R0,double w);
	PetscScalar inverse_r2_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots,double R0,double w);
	PetscScalar inverse_r_matrix_element(int i, int j, int degree, const std::vector<PetscScalar>& knots,double R0,double w);


	int save_bsplinee_basis(int n_basis, int degree, const std::vector<PetscScalar>& knots, int Nx, double xmax,double R0,double eta);
	PetscScalar R_x(PetscScalar x, double R0, double eta);
	std::vector<PetscScalar> R_k(const std::vector<PetscScalar>& knots, double R0, double eta);
	PetscScalar R_w(PetscScalar x, double w, double R0, double eta);

}
