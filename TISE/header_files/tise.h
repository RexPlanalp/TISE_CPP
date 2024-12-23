#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <petscmat.h>
#include <vector>
#include <string>
#include <complex>
namespace tise {
    PetscErrorCode construct_kinetic_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0,double eta);
    PetscErrorCode construct_inv_r2_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0, double eta);
    PetscErrorCode construct_overlap_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0, double eta);
    PetscErrorCode construct_inv_r_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0, double eta);
    PetscErrorCode solve_tise(int n_basis,int degree,const std::vector<std::complex<double>>& knots,int lmax, int nmax,double R0, double eta,double tolerance, PetscInt max_iter);
}

#endif // MATRIX_OPERATIONS_H
