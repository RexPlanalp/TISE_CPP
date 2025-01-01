#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <petscmat.h>
#include <vector>
#include <string>
#include <complex>
namespace tise {
    PetscErrorCode save_matrix(Mat A, const char *filename);
    PetscErrorCode construct_kinetic_matrix(Mat *A, PetscInt n_basis, PetscInt degree, const std::vector<PetscScalar>& knots,double R0,double eta);
    PetscErrorCode construct_inv_r2_matrix(Mat *A, PetscInt n_basis, PetscInt degree, const std::vector<PetscScalar>& knots,double R0, double eta);
    PetscErrorCode construct_overlap_matrix(Mat *A, PetscInt n_basis, PetscInt degree, const std::vector<PetscScalar>& knots,double R0, double eta);
    PetscErrorCode construct_inv_r_matrix(Mat *A, PetscInt n_basis, PetscInt degree, const std::vector<PetscScalar>& knots,double R0, double eta);
    PetscErrorCode solve_tise(PetscInt n_basis,PetscInt degree,const std::vector<PetscScalar>& knots,PetscInt lmax, PetscInt nmax,double R0, double eta,double tolerance, PetscInt  max_iter);
}

#endif // MATRIX_OPERATIONS_H
