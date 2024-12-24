#include <petscmat.h> // Include PETSc matrix header
#include <petscviewerhdf5.h>
#include <slepceps.h>
#include "basis.h"
#include <cstdlib>
#include <vector>
#include <fstream>
#include "tise.h"
#include <string>
#include <complex>
namespace tise 
{

    PetscErrorCode save_matrix(Mat A, const char *filename)
    {
    PetscErrorCode ierr;
    PetscViewer viewer;

    // Open a binary viewer in write mode
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

    // Write the matrix to the file in parallel
    ierr = MatView(A, viewer); CHKERRQ(ierr);

    // Clean up the viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    return ierr;
    }


    PetscErrorCode construct_kinetic_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0,double eta) 
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
                std::complex<double> matrix_element = basis::kinetic_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots,
                    R0,
                    eta
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

    PetscErrorCode construct_inv_r2_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0, double eta) 
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
                std::complex<double> matrix_element = basis::inverse_r2_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots,
                    R0,
                    eta
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
    
    PetscErrorCode construct_overlap_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0, double eta) 
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
                std::complex<double> matrix_element = basis::overlap_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots,
                    R0,
                    eta
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

    PetscErrorCode construct_inv_r_matrix(Mat *A, int n_basis, int degree, const std::vector<std::complex<double>>& knots,double R0, double eta) 
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
                std::complex<double> matrix_element = basis::inverse_r_matrix_element(
                    static_cast<int>(i), 
                    static_cast<int>(j), 
                    degree, 
                    knots,
                    R0,
                    eta
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

    PetscErrorCode solve_tise(int n_basis,int degree,const std::vector<std::complex<double>>& knots,int lmax, int nmax,double R0, double eta,double tolerance, PetscInt max_iter,bool EMBED)
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
        Mat H_total;
        Mat S_total;
        Mat I;
        Mat I_partial;
        Mat H_partial;

        ierr = construct_overlap_matrix(&S, n_basis, degree, knots,R0,eta); CHKERRQ(ierr);
        
        if (EMBED)
        {
            int total = 0;
            for (int l = 0; l<=lmax; l++)
            {
                total += 2*l + 1;
            }
            int total_size = total * n_basis;

            ierr = MatCreate(PETSC_COMM_WORLD, &I); CHKERRQ(ierr);
            ierr = MatSetSizes(I,PETSC_DECIDE,PETSC_DECIDE,total,total); CHKERRQ(ierr);
            ierr = MatSetFromOptions(I); CHKERRQ(ierr);
            PetscInt nnz_per_row_identity = 1;
            ierr = MatMPIAIJSetPreallocation(I,nnz_per_row_identity,NULL,nnz_per_row_identity,NULL); CHKERRQ(ierr);
            ierr = MatSetUp(I); CHKERRQ(ierr);
            PetscInt start_row,end_row;
            ierr = MatGetOwnershipRange(I,&start_row,&end_row); CHKERRQ(ierr);

            PetscInt nnz_per_row = 2*degree + 1;

            ierr = MatCreate(PETSC_COMM_WORLD, &H_total); CHKERRQ(ierr);
            ierr = MatSetSizes(H_total,PETSC_DECIDE,PETSC_DECIDE,total_size,total_size); CHKERRQ(ierr);
            ierr = MatSetFromOptions(H_total); CHKERRQ(ierr);
            ierr = MatMPIAIJSetPreallocation(H_total,nnz_per_row,NULL,nnz_per_row,NULL); CHKERRQ(ierr);
            ierr = MatSetUp(H_total); CHKERRQ(ierr);
            ierr = MatAssemblyBegin(H_total,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatAssemblyEnd(H_total,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

            ierr = MatCreate(PETSC_COMM_WORLD, &S_total); CHKERRQ(ierr);
            ierr = MatSetSizes(S_total,PETSC_DECIDE,PETSC_DECIDE,total_size,total_size); CHKERRQ(ierr);
            ierr = MatSetFromOptions(S_total); CHKERRQ(ierr);
            ierr = MatMPIAIJSetPreallocation(S_total,nnz_per_row,NULL,nnz_per_row,NULL); CHKERRQ(ierr);
            ierr = MatSetUp(S_total); CHKERRQ(ierr);
            ierr = MatAssemblyBegin(S_total,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatAssemblyEnd(S_total,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            
           
            for (PetscInt i = start_row; i< end_row; i++)
            {
                ierr = MatSetValue(I,i,i,1,INSERT_VALUES); CHKERRQ(ierr);
            }
           

            ierr = MatAssemblyBegin(I,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
            ierr = MatAssemblyEnd(I,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

            ierr = MatSeqAIJKron(I,S,MAT_INITIAL_MATRIX,&S_total); CHKERRQ(ierr);
            ierr = save_matrix(S_total,"S_total.dat"); CHKERRQ(ierr);

            ierr = MatDestroy(&S_total); CHKERRQ(ierr);
            ierr = MatDestroy(&I); CHKERRQ(ierr);
        }
        
       

        ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "tise_output.h5", FILE_MODE_WRITE, &viewTISE);

        ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
        ierr = EPSSetProblemType(eps, EPS_GNHEP); CHKERRQ(ierr);
        ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL); CHKERRQ(ierr);
        ierr = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERRQ(ierr);
        ierr = EPSSetTolerances(eps, tolerance, max_iter); CHKERRQ(ierr);
       


        ierr = construct_kinetic_matrix(&K, n_basis, degree, knots,R0,eta); CHKERRQ(ierr);
        ierr = construct_inv_r2_matrix(&Inv_r2, n_basis, degree,knots,R0,eta); CHKERRQ(ierr);
        ierr = construct_inv_r_matrix(&Inv_r,n_basis,degree,knots,R0,eta); CHKERRQ(ierr);

        for (int l=0; l<=lmax; ++l)
        {
            ierr = MatDuplicate(K, MAT_COPY_VALUES, &temp); CHKERRQ(ierr);
            ierr = MatAXPY(temp, l*(l+1)*0.5,Inv_r2,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
            ierr = MatAXPY(temp,-1.0,Inv_r,SAME_NONZERO_PATTERN); CHKERRQ(ierr);

            if (EMBED)
            {   

                
                int total = 0;
                for (int l = 0; l<=lmax; l++)
                {
                    total += 2*l + 1;
                }
                int total_size = total * n_basis;


                int start = 0;
                for (int l_temp = 0; l_temp<l; l_temp++)
                {
                    start += 2*l_temp + 1;
                }
                int end = start + 2*l+1;


                PetscInt nnz_per_row = 2*degree +1;
                ierr = MatCreate(PETSC_COMM_WORLD, &H_partial); CHKERRQ(ierr);
                ierr = MatSetSizes(H_partial,PETSC_DECIDE,PETSC_DECIDE,total_size,total_size); CHKERRQ(ierr);
                ierr = MatSetFromOptions(H_total); CHKERRQ(ierr);
                ierr = MatMPIAIJSetPreallocation(H_partial,nnz_per_row,NULL,nnz_per_row,NULL); CHKERRQ(ierr);
                ierr = MatSetUp(H_partial); CHKERRQ(ierr);
                ierr = MatAssemblyBegin(H_partial,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
                ierr = MatAssemblyEnd(H_partial,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

                ierr = MatCreate(PETSC_COMM_WORLD, &I_partial); CHKERRQ(ierr);
                ierr = MatSetSizes(I_partial,PETSC_DECIDE,PETSC_DECIDE,total,total); CHKERRQ(ierr);
                ierr = MatSetFromOptions(I_partial); CHKERRQ(ierr);
                PetscInt nnz_per_row_identity = 1;
                ierr = MatMPIAIJSetPreallocation(I_partial,nnz_per_row_identity,NULL,nnz_per_row_identity,NULL); CHKERRQ(ierr);
                ierr = MatSetUp(I_partial); CHKERRQ(ierr);
                PetscInt start_row,end_row;
                ierr = MatGetOwnershipRange(I_partial,&start_row,&end_row); CHKERRQ(ierr);

                for (PetscInt i = start_row; i< end_row; i++)
                {
                    if ((i >= start && i < end))
                    {
                        ierr = MatSetValue(I_partial,i,i,1,INSERT_VALUES); CHKERRQ(ierr);
                    }
                    else
                    {
                        continue;
                    }
                }

                ierr = MatAssemblyBegin(I_partial,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
                ierr = MatAssemblyEnd(I_partial,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

                ierr = MatSeqAIJKron(I_partial,temp,MAT_INITIAL_MATRIX,&H_partial); CHKERRQ(ierr);

                ierr = MatAXPY(H_total,1.0,H_partial,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
                ierr = MatDestroy(&H_partial); CHKERRQ(ierr);
                ierr = MatDestroy(&I_partial); CHKERRQ(ierr);
            }

            int num_of_energies = nmax - l;
            if (num_of_energies <= 0)
            {
                continue;
            }


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
        if (EMBED)
        {
            ierr = save_matrix(H_total,"H_total.dat"); CHKERRQ(ierr);
            ierr = MatDestroy(&H_total); CHKERRQ(ierr);
        }
        


        return 0;

    }
} // namespace tise
