#include <petscmat.h> // Include PETSc matrix header
#include <petscviewerhdf5.h>
#include <slepceps.h>
#include "basis.h"
#include <cstdlib>
#include <vector>
#include <fstream>
#include "tise.h"
#include <string>
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
    
    PetscErrorCode construct_inv_r_matrix(Mat *A, int n_basis, int degree, const std::vector<double>& knots) 
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
                double matrix_element = basis::inverse_r_matrix_element(
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

    PetscErrorCode solve_tise(int n_basis,int degree,const std::vector<double>& knots,int lmax, int nmax)
    {
        PetscErrorCode ierr;
        PetscViewer viewTISE;
        Mat K;
        Mat Inv_r2;
        Mat Inv_r;
        Mat S;
        Mat temp;
        EPS eps;
        PetscInt nconv;
        

        double tolerance = 1e-15;
        PetscInt max_iter = 5000;


        ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "tise_output.h5", FILE_MODE_WRITE, &viewTISE);

        ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
        ierr = EPSSetProblemType(eps, EPS_GNHEP); CHKERRQ(ierr);
        ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL); CHKERRQ(ierr);
        ierr = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERRQ(ierr);
        ierr = EPSSetTolerances(eps, tolerance, max_iter); CHKERRQ(ierr);
       


        ierr = construct_kinetic_matrix(&K, n_basis, degree, knots); CHKERRQ(ierr);
        ierr = construct_overlap_matrix(&S, n_basis, degree, knots); CHKERRQ(ierr);
        ierr = construct_inv_r2_matrix(&Inv_r2, n_basis, degree,knots); CHKERRQ(ierr);
        ierr = construct_inv_r_matrix(&Inv_r,n_basis,degree,knots); CHKERRQ(ierr);

        for (int l=0; l<=lmax; ++l)
        {
            int num_of_energies = nmax - l;
            if (num_of_energies <= 0)
            {
                continue;
            }

            ierr = MatDuplicate(K, MAT_COPY_VALUES, &temp); CHKERRQ(ierr);
            ierr = MatAXPY(temp, l*(l+1)*0.5,Inv_r2,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
            ierr = MatAXPY(temp,-1.0,Inv_r,SAME_NONZERO_PATTERN); CHKERRQ(ierr);

            ierr = EPSSetOperators(eps, temp, S); CHKERRQ(ierr);
            ierr = EPSSetDimensions(eps,num_of_energies,PETSC_DEFAULT,PETSC_DEFAULT);
            ierr = EPSSolve(eps); CHKERRQ(ierr);
            ierr = EPSGetConverged(eps, &nconv); CHKERRQ(ierr);
            PetscPrintf(PETSC_COMM_WORLD, "Eigenvalues Requested: %d, Eigenvalues Converged: %d\n", num_of_energies, nconv);


            for (PetscInt i = 0; i < nconv; ++i)
            {
                PetscScalar eigenvalue;
                ierr = EPSGetEigenvalue(eps, i, &eigenvalue, NULL); CHKERRQ(ierr);

                PetscReal real_part;
                real_part = PetscRealPart(eigenvalue);

                if (real_part>0)
                {
                    continue;
                }
                
                // Retrieve eigenvalue
                ierr = EPSGetEigenvalue(eps, i, &eigenvalue, NULL); CHKERRQ(ierr);

                // Save eigenvalue as a scalar
                std::string eigenvalue_name = std::string("E_") + std::to_string(i + l + 1) + "_" + std::to_string(l);
                ierr = PetscViewerHDF5PushGroup(viewTISE, "/eigenvalues"); CHKERRQ(ierr);

                Vec eigenvalue_vec;
                ierr = VecCreate(PETSC_COMM_WORLD, &eigenvalue_vec); CHKERRQ(ierr);
                ierr = VecSetSizes(eigenvalue_vec, PETSC_DECIDE, 1); CHKERRQ(ierr);
                ierr = VecSetFromOptions(eigenvalue_vec); CHKERRQ(ierr);
                ierr = VecSetValue(eigenvalue_vec, 0, eigenvalue, INSERT_VALUES); CHKERRQ(ierr);
                ierr = VecAssemblyBegin(eigenvalue_vec); CHKERRQ(ierr);
                ierr = VecAssemblyEnd(eigenvalue_vec); CHKERRQ(ierr);

                ierr = PetscObjectSetName((PetscObject)eigenvalue_vec, eigenvalue_name.c_str()); CHKERRQ(ierr);
                ierr = VecView(eigenvalue_vec, viewTISE); CHKERRQ(ierr);
                ierr = PetscViewerHDF5PopGroup(viewTISE); CHKERRQ(ierr);
                ierr = VecDestroy(&eigenvalue_vec); CHKERRQ(ierr);

                // Retrieve and save eigenvector
                Vec eigenvector;
                ierr = MatCreateVecs(temp, &eigenvector, NULL); CHKERRQ(ierr);
                ierr = EPSGetEigenvector(eps, i, eigenvector, NULL); CHKERRQ(ierr);

                std::string eigenvector_name = std::string("psi_l_") + std::to_string(i + l + 1) + "_" + std::to_string(l);
                ierr = PetscViewerHDF5PushGroup(viewTISE, "/eigenvectors"); CHKERRQ(ierr);
                ierr = PetscObjectSetName((PetscObject)eigenvector, eigenvector_name.c_str()); CHKERRQ(ierr);
                ierr = VecView(eigenvector, viewTISE); CHKERRQ(ierr);
                ierr = PetscViewerHDF5PopGroup(viewTISE); CHKERRQ(ierr);

                // Cleanup
                ierr = VecDestroy(&eigenvector); CHKERRQ(ierr);
            }
        
        }

        ierr = PetscViewerDestroy(&viewTISE); CHKERRQ(ierr);
        ierr = EPSDestroy(&eps); CHKERRQ(ierr);
        ierr = MatDestroy(&K); CHKERRQ(ierr);
        ierr = MatDestroy(&Inv_r2); CHKERRQ(ierr);
        ierr = MatDestroy(&Inv_r); CHKERRQ(ierr);
        ierr = MatDestroy(&S); CHKERRQ(ierr);
        ierr = MatDestroy(&temp); CHKERRQ(ierr);


        return 0;

    }



    

} // namespace tise
