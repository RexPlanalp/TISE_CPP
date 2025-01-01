#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <functional>
#include <petscsys.h>

namespace basis 
{	
	std::vector<PetscScalar> linear_knots(PetscInt n_basis, PetscInt degree, double rmax);
	PetscScalar B(PetscInt i, PetscInt degree, const std::vector<PetscScalar>& knots, PetscScalar x);
	PetscScalar dB(PetscInt i, PetscInt degree, const std::vector<PetscScalar>& knots, PetscScalar x);


	PetscScalar compute_matrix_element(PetscInt i, PetscInt j, PetscInt degree, const std::vector<PetscScalar>& knots, std::function<PetscScalar(int, int, const std::vector<PetscScalar>&, PetscScalar)> integrand, double R0, double eta);
	PetscScalar overlap_matrix_element(PetscInt i, PetscInt j, PetscInt degree, const std::vector<PetscScalar>& knots, double R0, double eta);
	PetscScalar kinetic_matrix_element(PetscInt i, PetscInt j, PetscInt degree, const std::vector<PetscScalar>& knots,double R0,double eta);
	PetscScalar inverse_r2_matrix_element(PetscInt i, PetscInt j, PetscInt degree, const std::vector<PetscScalar>& knots,double R0,double eta);
	PetscScalar inverse_r_matrix_element(PetscInt i, PetscInt j, PetscInt degree, const std::vector<PetscScalar>& knots,double R0,double eta);


	PetscErrorCode save_bsplinee_basis(PetscInt n_basis, PetscInt degree, const std::vector<PetscScalar>& knots, PetscInt Nx, double xmax,double R0,double eta);
	PetscScalar R_x(PetscScalar x, double R0, double eta);
	std::vector<PetscScalar> R_k(const std::vector<PetscScalar>& knots, double R0, double eta);
	PetscScalar R_w(PetscScalar x, double w, double R0, double eta);

}
