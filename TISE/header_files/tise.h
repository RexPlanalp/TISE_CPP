#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <petscmat.h>
#include <vector>

namespace tise {
    PetscErrorCode construct_kinetic_matrix(Mat *A, int n_basis, int degree,const std::vector<double>& knots);
    PetscErrorCode construct_inv_r2_matrix(Mat *A, int n_basis, int degree, const std::vector<double>& knots);
    PetscErrorCode construct_overlap_matrix(Mat *A, int n_basis, int degree, const std::vector<double>& knots);
}

#endif // MATRIX_OPERATIONS_H
