#include <petscmat.h> // Include PETSc matrix header
#include "basis.h"
#include <cstdlib>
#include <vector>

namespace tise 
{

    PetscErrorCode construct_kinetic_matrix(Mat *A, int n_basis, int degree, const std::vector<double>& knots) 
    {
        PetscInt p_n_basis = static_cast<PetscInt>(n_basis);
        PetscInt p_degree = static_cast<PetscInt>(degree);

        PetscErrorCode ierr;

        // Create the matrix
        ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
        ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, p_n_basis, p_n_basis); CHKERRQ(ierr);
        ierr = MatSetFromOptions(*A); CHKERRQ(ierr);

        // Preallocate memory for nonzero entries
        PetscInt nnz_per_row = 2 * p_degree + 1; // Number of nonzeros per row
        ierr = MatMPIAIJSetPreallocation(*A, nnz_per_row, NULL, nnz_per_row, NULL); CHKERRQ(ierr);

        // Set up the matrix
        ierr = MatSetUp(*A); CHKERRQ(ierr);

        // Get the range of rows owned by the current process
        PetscInt start_row, end_row;
        ierr = MatGetOwnershipRange(*A, &start_row, &end_row); CHKERRQ(ierr);

        // Precompute the degree-based band width
        PetscInt band_width = p_degree + 1;

        // Iterate over locally owned rows
        for (PetscInt i = start_row; i < end_row; i++) 
        {
            // Set values only within the band for current row
            PetscInt col_start = std::max(static_cast<PetscInt>(0), i - band_width + 1);
            PetscInt col_end = std::min(p_n_basis, i + band_width); // Exclusive

            for (PetscInt j = col_start; j < col_end; j++) 
            {
                // Compute the matrix element
                double matrix_element = basis::kinetic_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots
                );

                // Set the value in the matrix
                ierr = MatSetValue(*A, i, j, matrix_element, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
        
        // Assemble the matrix
        ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        return 0; // Return success
    }

    PetscErrorCode construct_inv_r2_matrix(Mat *A, int n_basis, int degree, const std::vector<double>& knots) 
    {
        PetscInt p_n_basis = static_cast<PetscInt>(n_basis);
        PetscInt p_degree = static_cast<PetscInt>(degree);

        PetscErrorCode ierr;

        // Create the matrix
        ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
        ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, p_n_basis, p_n_basis); CHKERRQ(ierr);
        ierr = MatSetFromOptions(*A); CHKERRQ(ierr);

        // Preallocate memory for nonzero entries
        PetscInt nnz_per_row = 2 * p_degree + 1; // Number of nonzeros per row
        ierr = MatMPIAIJSetPreallocation(*A, nnz_per_row, NULL, nnz_per_row, NULL); CHKERRQ(ierr);

        // Set up the matrix
        ierr = MatSetUp(*A); CHKERRQ(ierr);

        // Get the range of rows owned by the current process
        PetscInt start_row, end_row;
        ierr = MatGetOwnershipRange(*A, &start_row, &end_row); CHKERRQ(ierr);

        // Precompute the degree-based band width
        PetscInt band_width = p_degree + 1;

        // Iterate over locally owned rows
        for (PetscInt i = start_row; i < end_row; i++) 
        {
            // Set values only within the band for current row
            PetscInt col_start = std::max(static_cast<PetscInt>(0), i - band_width + 1);
            PetscInt col_end = std::min(p_n_basis, i + band_width); // Exclusive

            for (PetscInt j = col_start; j < col_end; j++) 
            {
                // Compute the matrix element
                double matrix_element = basis::inverse_r2_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots
                );

                // Set the value in the matrix
                ierr = MatSetValue(*A, i, j, matrix_element, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
        
        // Assemble the matrix
        ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        return 0; // Return success
    }
    
    PetscErrorCode construct_overlap_matrix(Mat *A, int n_basis, int degree, const std::vector<double>& knots) 
    {
        PetscInt p_n_basis = static_cast<PetscInt>(n_basis);
        PetscInt p_degree = static_cast<PetscInt>(degree);

        PetscErrorCode ierr;

        // Create the matrix
        ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
        ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, p_n_basis, p_n_basis); CHKERRQ(ierr);
        ierr = MatSetFromOptions(*A); CHKERRQ(ierr);

        // Preallocate memory for nonzero entries
        PetscInt nnz_per_row = 2 * p_degree + 1; // Number of nonzeros per row
        ierr = MatMPIAIJSetPreallocation(*A, nnz_per_row, NULL, nnz_per_row, NULL); CHKERRQ(ierr);

        // Set up the matrix
        ierr = MatSetUp(*A); CHKERRQ(ierr);

        // Get the range of rows owned by the current process
        PetscInt start_row, end_row;
        ierr = MatGetOwnershipRange(*A, &start_row, &end_row); CHKERRQ(ierr);

        // Precompute the degree-based band width
        PetscInt band_width = p_degree + 1;

        // Iterate over locally owned rows
        for (PetscInt i = start_row; i < end_row; i++) 
        {
            // Set values only within the band for current row
            PetscInt col_start = std::max(static_cast<PetscInt>(0), i - band_width + 1);
            PetscInt col_end = std::min(p_n_basis, i + band_width); // Exclusive

            for (PetscInt j = col_start; j < col_end; j++) 
            {
                // Compute the matrix element
                double matrix_element = basis::overlap_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots
                );

                // Set the value in the matrix
                ierr = MatSetValue(*A, i, j, matrix_element, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
        
        // Assemble the matrix
        ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        return 0; // Return success
    }



} // namespace tise
