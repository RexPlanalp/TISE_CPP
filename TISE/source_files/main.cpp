// #include "tise.h"
// #include <petscsys.h>
// #include <slepceps.h>
// #include <vector>
// #include "basis.h"
// #include <mpi.h>
// #include <complex>

// PetscInt main(PetscInt argc, char **argv) {
    
//     PetscErrorCode ierr;
//     ierr = SlepcInitialize(&argc, &argv, NULL,NULL); CHKERRQ(ierr);

//     PetscInt n_basis = 3000;
//     PetscInt degree = 6;
//     double rmax = 3000;
//     PetscInt lmax = 75;
//     PetscInt nmax = 10;
//     double R0 = 3000;
//     double eta=  0;
//     double tise_tolerance = 1E-15;
//     PetscInt  tise_mat_iter = 3000;
//     bool SAVE_BSPLINES = false;

//     std::vector<PetscScalar> knots = basis::linear_knots(n_basis, degree, rmax);

//     if (SAVE_BSPLINES)
//     {
//         PetscInt Nx = rmax * 10;
//         ierr = basis::save_bsplinee_basis(n_basis,degree,knots,Nx,rmax,R0,eta); CHKERRQ(ierr);
//     }
  
//     PetscPrintf(PETSC_COMM_WORLD, "Solving TISE\n");
//     double start_time = MPI_Wtime();
//     tise::solve_tise(n_basis,degree,knots,lmax,nmax,R0,eta,tise_tolerance,tise_mat_iter);
//     double end_time = MPI_Wtime();
//     double elapsed_time = end_time - start_time;
//     PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

//     ierr = SlepcFinalize(); CHKERRQ(ierr);
//     return 0;

// }

// #include <petscksp.h>
// #include <omp.h> // OpenMP for parallelism

// int main(int argc, char **argv) {
//     PetscErrorCode ierr;
//     ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

//     const PetscInt nA = 6000, nB = 6000; // Matrix dimensions
//     const PetscInt n = nA * nB;

//     Mat A, B;
//     ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nA, nA, 4, NULL, &A); CHKERRQ(ierr);
//     ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nB, nB, 3, NULL, &B); CHKERRQ(ierr);

//     // Fill matrices A and B
//     for (PetscInt i = 0; i < nA; ++i) {
//         for (PetscInt j = 0; j < nA; ++j) {
//             if (j == i || j == (i + 1) % nA || j == (i + 2) % nA || j == (i + 3) % nA) {
//                 ierr = MatSetValue(A, i, j, (i + 1) * (j + 1), INSERT_VALUES); CHKERRQ(ierr);
//             }
//         }
//     }
//     for (PetscInt i = 0; i < nB; ++i) {
//         for (PetscInt j = 0; j < nB; ++j) {
//             if (j == i || j == (i + 1) % nB || j == (i + 2) % nB) {
//                 ierr = MatSetValue(B, i, j, (i + 1) + (j + 1), INSERT_VALUES); CHKERRQ(ierr);
//             }
//         }
//     }
//     ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//     ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//     ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//     ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

//     Vec x, y1, y2;
//     ierr = VecCreateSeq(PETSC_COMM_SELF, n, &x); CHKERRQ(ierr);
//     ierr = VecDuplicate(x, &y1); CHKERRQ(ierr);
//     ierr = VecDuplicate(x, &y2); CHKERRQ(ierr);
//     for (PetscInt i = 0; i < n; ++i) {
//         ierr = VecSetValue(x, i, i + 1.0, INSERT_VALUES); CHKERRQ(ierr);
//     }
//     ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
//     ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

//     PetscLogDouble t_start, t_end;

//     // Full tensor product multiplication
//     Mat AB;
//     ierr = MatSeqAIJKron(A, B, MAT_INITIAL_MATRIX, &AB); CHKERRQ(ierr);
//     ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
//     t_start = MPI_Wtime();
//     ierr = MatMult(AB, x, y1); CHKERRQ(ierr);
//     t_end = MPI_Wtime();
//     PetscPrintf(PETSC_COMM_SELF, "Time for full matrix multiplication: %f seconds\n", t_end - t_start);

//     PetscScalar sum_y1;
//     ierr = VecSum(y1, &sum_y1); CHKERRQ(ierr);

//     // Optimized 3-loop implementation
//     PetscInt nrowsA, nrowsB;
//     const PetscInt *iA, *jA, *iB, *jB;
//     PetscBool done;

//     // Retrieve CSR data for matrix A
//     ierr = MatGetRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &nrowsA, &iA, &jA, &done); CHKERRQ(ierr);
//     if (!done) {
//         PetscPrintf(PETSC_COMM_SELF, "MatGetRowIJ did not complete successfully for matrix A.\n");
//         PetscFinalize();
//         return 1;
//     }

//     // Retrieve CSR data for matrix B
//     ierr = MatGetRowIJ(B, 0, PETSC_FALSE, PETSC_FALSE, &nrowsB, &iB, &jB, &done); CHKERRQ(ierr);
//     if (!done) {
//         PetscPrintf(PETSC_COMM_SELF, "MatGetRowIJ did not complete successfully for matrix B.\n");
//         PetscFinalize();
//         return 1;
//     }

//     const PetscScalar *x_array;
//     PetscScalar *y2_array;
//     ierr = VecGetArrayRead(x, &x_array); CHKERRQ(ierr);
//     ierr = VecGetArray(y2, &y2_array); CHKERRQ(ierr);

//     PetscArrayzero(y2_array, n);

//     ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
//     t_start = MPI_Wtime();
//     #pragma omp parallel for
//     for (PetscInt i = 0; i < nA; ++i) {
//         PetscInt rowOffset = i * nB;
//         PetscInt startA = iA[i], endA = iA[i + 1];
//         for (PetscInt kA = startA; kA < endA; ++kA) {
//             PetscInt colA = jA[kA];
//             PetscScalar vA = (PetscScalar)((i + 1) * (colA + 1));
//             PetscInt colOffset = colA * nB;
//             for (PetscInt j = 0; j < nB; ++j) {
//                 PetscInt startB = iB[j], endB = iB[j + 1];
//                 PetscScalar partialSum = 0.0;
//                 for (PetscInt kB = startB; kB < endB; ++kB) {
//                     PetscInt colB = jB[kB];
//                     partialSum += ((j + 1) + (colB + 1)) * x_array[colOffset + colB];
//                 }
//                 y2_array[rowOffset + j] += vA * partialSum;
//             }
//         }
//     }
//     ierr = VecRestoreArrayRead(x, &x_array); CHKERRQ(ierr);
//     ierr = VecRestoreArray(y2, &y2_array); CHKERRQ(ierr);

//     ierr = MatRestoreRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &nrowsA, &iA, &jA, &done); CHKERRQ(ierr);
//     ierr = MatRestoreRowIJ(B, 0, PETSC_FALSE, PETSC_FALSE, &nrowsB, &iB, &jB, &done); CHKERRQ(ierr);

//     t_end = MPI_Wtime();
//     PetscPrintf(PETSC_COMM_SELF, "Time for optimized computation using CSR (3 loops): %f seconds\n", t_end - t_start);

//     PetscScalar sum_y2;
//     ierr = VecSum(y2, &sum_y2); CHKERRQ(ierr);

//     PetscReal rel_error = PetscAbsReal(sum_y1 - sum_y2) / PetscAbsReal(sum_y1);
//     if (rel_error < 1e-12) {
//         PetscPrintf(PETSC_COMM_SELF, "The results match within tolerance!\n");
//     } else {
//         PetscPrintf(PETSC_COMM_SELF, "The results do NOT match (relative error: %e).\n", rel_error);
//     }

//     ierr = MatDestroy(&A); CHKERRQ(ierr);
//     ierr = MatDestroy(&B); CHKERRQ(ierr);
//     ierr = MatDestroy(&AB); CHKERRQ(ierr);
//     ierr = VecDestroy(&x); CHKERRQ(ierr);
//     ierr = VecDestroy(&y1); CHKERRQ(ierr);
//     ierr = VecDestroy(&y2); CHKERRQ(ierr);

//     ierr = PetscFinalize(); CHKERRQ(ierr);
//     return 0;
// }





// #include <petscksp.h>
// #include <vector>
// #include <iostream>
// #include <unordered_map>

// // Function to compute y = (A kron B) * x
// void kronSpMV(const Mat &A, const Mat &B, const Vec &x, Vec &y) {
//     PetscInt N, M;
//     MatGetSize(A, &N, nullptr); // Get dimension of A (N x N)
//     MatGetSize(B, &M, nullptr); // Get dimension of B (M x M)
//     PetscInt xSize, ySize;
//     VecGetSize(x, &xSize);
//     VecGetSize(y, &ySize);

//     if (xSize != N * M || ySize != N * M) {
//         std::cerr << "Error: Mismatch in dimensions of vector x or y!" << std::endl;
//         return;
//     }

//     PetscInt local_start, local_end;
//     VecGetOwnershipRange(y, &local_start, &local_end);

//     const PetscScalar *xArray;
//     VecGetArrayRead(x, &xArray);

//     for (PetscInt i = 0; i < N; ++i) {
//         PetscInt AStart, AEnd;
//         const PetscInt *ACols;
//         const PetscScalar *AValues;
//         MatGetRow(A, i, &AEnd, &ACols, &AValues);

//         for (PetscInt k = 0; k < AEnd; ++k) {
//             PetscInt p = ACols[k];         // Column index of A
//             PetscScalar valA = AValues[k]; // Value of A[i, p]

//             for (PetscInt j = 0; j < M; ++j) {
//                 PetscScalar sum_b = 0.0;

//                 PetscInt BStart, BEnd;
//                 const PetscInt *BCols;
//                 const PetscScalar *BValues;
//                 MatGetRow(B, j, &BEnd, &BCols, &BValues);

//                 // Compute contributions to sum_b
//                 for (PetscInt l = 0; l < BEnd; ++l) {
//                     PetscInt q = BCols[l];
//                     PetscScalar valB = BValues[l];
//                     PetscInt globalIndex = p * M + q;

//                     PetscScalar value = 0.0;
//                     if (globalIndex >= local_start && globalIndex < local_end) {
//                         // Use locally owned part of the vector
//                         value = xArray[globalIndex - local_start];
//                     } else {
//                         // For non-local indices, add scatter logic here if needed
//                         // or assume they were pre-gathered using VecScatter.
//                         continue;
//                     }

//                     sum_b += valB * value;
//                 }

//                 MatRestoreRow(B, j, &BEnd, &BCols, &BValues);

//                 // Update global y using VecSetValue
//                 PetscInt globalIndex = i * M + j;
//                 VecSetValue(y, globalIndex, valA * sum_b, ADD_VALUES);
//             }
//         }
//         MatRestoreRow(A, i, &AEnd, &ACols, &AValues);
//     }

//     VecRestoreArrayRead(x, &xArray);

//     // Assemble the global vector y
//     VecAssemblyBegin(y);
//     VecAssemblyEnd(y);
// }

// int main(int argc, char **args) {
//     PetscInitialize(&argc, &args, nullptr, nullptr);

//     // Define small matrix A (2x2) and B (3x3)
//     Mat A, B;
//     Vec x, y;

//     PetscInt nA = 2, nB = 3; // Dimensions
//     MatCreateSeqAIJ(PETSC_COMM_SELF, nA, nA, 2, nullptr, &A);
//     MatCreateSeqAIJ(PETSC_COMM_SELF, nB, nB, 3, nullptr, &B);

//     // Fill A
//     MatSetValue(A, 0, 0, 2.0, INSERT_VALUES);
//     MatSetValue(A, 0, 1, 1.0, INSERT_VALUES);
//     MatSetValue(A, 1, 0, 3.0, INSERT_VALUES);
//     MatSetValue(A, 1, 1, 4.0, INSERT_VALUES);
//     MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

//     // Fill B
//     MatSetValue(B, 0, 0, 1.0, INSERT_VALUES);
//     MatSetValue(B, 0, 1, 0.5, INSERT_VALUES);
//     MatSetValue(B, 0, 2, 0.2, INSERT_VALUES);
//     MatSetValue(B, 1, 0, 0.1, INSERT_VALUES);
//     MatSetValue(B, 1, 1, 1.5, INSERT_VALUES);
//     MatSetValue(B, 1, 2, 0.3, INSERT_VALUES);
//     MatSetValue(B, 2, 0, 0.4, INSERT_VALUES);
//     MatSetValue(B, 2, 1, 0.2, INSERT_VALUES);
//     MatSetValue(B, 2, 2, 2.0, INSERT_VALUES);
//     MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

//     // Initialize vector x with random values
//     PetscInt size = nA * nB;
//     VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &x);
//     VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &y);

//     PetscMPIInt rank;
//     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
//     PetscInt local_start, local_end;
//     VecGetOwnershipRange(x, &local_start, &local_end);

//     for (PetscInt i = local_start; i < local_end; ++i) {
//         PetscScalar val = static_cast<PetscScalar>(rand()) / RAND_MAX;
//         VecSetValue(x, i, val, INSERT_VALUES);
//     }
//     VecAssemblyBegin(x);
//     VecAssemblyEnd(x);

//     double start = MPI_Wtime();
//     // Compute y = (A kron B) * x
//     kronSpMV(A, B, x, y);
//     double end = MPI_Wtime();
//     double elapsed = end-start;
//     std::cout<<elapsed<<std::endl;


//     // Print y
//     PetscViewer viewer;
//     PetscViewerASCIIGetStdout(PETSC_COMM_WORLD, &viewer);
//     PetscViewerASCIIPrintf(viewer, "Process [%d]\n", rank);
//     VecView(y, viewer);

//     // Clean up
//     MatDestroy(&A);
//     MatDestroy(&B);
//     VecDestroy(&x);
//     VecDestroy(&y);

//     PetscFinalize();
//     return 0;
// }



#include <petscksp.h>
#include <vector>
#include <unordered_map>
#include <iostream>

// WORKING?
void kronSpMV(const Mat &A, const Mat &B, const Vec &x, Vec &y) {
    PetscInt N, M;
    MatGetSize(A, &N, nullptr); // Get dimension of A (N x N)
    MatGetSize(B, &M, nullptr); // Get dimension of B (M x M)
    PetscInt local_start, local_end;
    VecGetOwnershipRange(y, &local_start, &local_end);

    const PetscScalar *xArray;
    PetscScalar *yArray; // Declare yArray
    VecGetArrayRead(x, &xArray);
    VecGetArray(y, &yArray);

    for (PetscInt i = 0; i < N; ++i) {
        PetscInt AStart, AEnd;
        const PetscInt *ACols;
        const PetscScalar *AValues;
        MatGetRow(A, i, &AEnd, &ACols, &AValues);

        for (PetscInt k = 0; k < AEnd; ++k) {
            PetscInt p = ACols[k];         // Column index of A
            PetscScalar valA = AValues[k]; // Value of A[i, p]

            for (PetscInt j = 0; j < M; ++j) {
                if (i * M + j < local_start || i * M + j >= local_end) {
                    continue; // Skip non-local parts of y
                }

                PetscScalar sum_b = 0.0;

                PetscInt BStart, BEnd;
                const PetscInt *BCols;
                const PetscScalar *BValues;
                MatGetRow(B, j, &BEnd, &BCols, &BValues);

                std::vector<PetscInt> neededIndices;
                std::unordered_map<PetscInt, PetscInt> globalToLocalMap;
                PetscInt localCounter = 0;

                // Identify non-local indices of x
                for (PetscInt l = 0; l < BEnd; ++l) {
                    PetscInt q = BCols[l];
                    PetscInt globalIndex = p * M + q;

                    if (globalIndex < local_start || globalIndex >= local_end) {
                        neededIndices.push_back(globalIndex);
                        globalToLocalMap[globalIndex] = localCounter++;
                    }
                }

                // Scatter non-local parts of x
                Vec localVec = nullptr;
                const PetscScalar *localArray = nullptr;
                if (!neededIndices.empty()) {
                    IS indexSet;
                    ISCreateGeneral(PETSC_COMM_SELF, neededIndices.size(), neededIndices.data(), PETSC_COPY_VALUES, &indexSet);

                    VecCreateSeq(PETSC_COMM_SELF, neededIndices.size(), &localVec);

                    VecScatter scatter;
                    VecScatterCreate(x, indexSet, localVec, nullptr, &scatter);
                    VecScatterBegin(scatter, x, localVec, INSERT_VALUES, SCATTER_FORWARD);
                    VecScatterEnd(scatter, x, localVec, INSERT_VALUES, SCATTER_FORWARD);

                    VecGetArrayRead(localVec, &localArray);

                    ISDestroy(&indexSet);
                    VecScatterDestroy(&scatter);
                }

                // Compute sum_b using local and scattered parts of x
                for (PetscInt l = 0; l < BEnd; ++l) {
                    PetscInt q = BCols[l];
                    PetscInt globalIndex = p * M + q;

                    PetscScalar value = 0.0;
                    if (globalIndex >= local_start && globalIndex < local_end) {
                        // Use local part of x
                        value = xArray[globalIndex - local_start];
                    } else {
                        // Use scattered part of x
                        value = localArray[globalToLocalMap[globalIndex]];
                    }

                    PetscScalar valB = BValues[l];
                    sum_b += valB * value;
                }

                if (localVec) {
                    VecRestoreArrayRead(localVec, &localArray);
                    VecDestroy(&localVec);
                }

                MatRestoreRow(B, j, &BEnd, &BCols, &BValues);

                // Update local part of y
                yArray[i * M + j - local_start] += valA * sum_b;
            }
        }
        MatRestoreRow(A, i, &AEnd, &ACols, &AValues);
    }

    VecRestoreArrayRead(x, &xArray);
    VecRestoreArray(y, &yArray); // Restore yArray
}





int main(int argc, char **args) {
    PetscInitialize(&argc, &args, nullptr, nullptr);

    // Define small matrix A (2x2) and B (3x3)
    Mat A, B;
    Vec x, y;

    PetscInt nA = 2, nB = 3; // Dimensions
    MatCreateSeqAIJ(PETSC_COMM_SELF, nA, nA, 2, nullptr, &A);
    MatCreateSeqAIJ(PETSC_COMM_SELF, nB, nB, 3, nullptr, &B);

    // Fill A
    MatSetValue(A, 0, 0, 2.0, INSERT_VALUES);
    MatSetValue(A, 0, 1, 1.0, INSERT_VALUES);
    MatSetValue(A, 1, 0, 3.0, INSERT_VALUES);
    MatSetValue(A, 1, 1, 4.0, INSERT_VALUES);
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    // Fill B
    MatSetValue(B, 0, 0, 1.0, INSERT_VALUES);
    MatSetValue(B, 0, 1, 0.5, INSERT_VALUES);
    MatSetValue(B, 0, 2, 0.2, INSERT_VALUES);
    MatSetValue(B, 1, 0, 0.1, INSERT_VALUES);
    MatSetValue(B, 1, 1, 1.5, INSERT_VALUES);
    MatSetValue(B, 1, 2, 0.3, INSERT_VALUES);
    MatSetValue(B, 2, 0, 0.4, INSERT_VALUES);
    MatSetValue(B, 2, 1, 0.2, INSERT_VALUES);
    MatSetValue(B, 2, 2, 2.0, INSERT_VALUES);
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);

    // Initialize vector x with random values
    PetscInt size = nA * nB;
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &x);
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &y);

    PetscInt local_start, local_end;
    VecGetOwnershipRange(x, &local_start, &local_end);

    for (PetscInt i = local_start; i < local_end; ++i) {
        PetscScalar val = static_cast<PetscScalar>(i + 1); // Deterministic initialization
        VecSetValue(x, i, val, INSERT_VALUES);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);


    double start = MPI_Wtime();
    // Compute y = (A kron B) * x
    kronSpMV(A, B, x, y);
    double end = MPI_Wtime();
    double elapsed = end-start;
    std::cout << elapsed << std::endl;
    PetscPrintf(PETSC_COMM_WORLD, "Time in ms %.3f \n", elapsed*1000);

    // Print y
    VecView(y, PETSC_VIEWER_STDOUT_WORLD);

    // Clean up
    MatDestroy(&A);
    MatDestroy(&B);
    VecDestroy(&x);
    VecDestroy(&y);

    PetscFinalize();
    return 0;

    // // Create the Kronecker product matrix
    // Mat K;
    // MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, nB * nA, nullptr, &K);
    // MatSeqAIJKron(A,B,MAT_INITIAL_MATRIX,&K);

    // // Measure time for K * x
    // double start = MPI_Wtime();
    // MatMult(K, x, y);
    // double end = MPI_Wtime();
    // double elapsed = end - start;

    // // Print timing in ms
    // PetscPrintf(PETSC_COMM_SELF, "Time to multiply K * x: %.3f ms\n", elapsed * 1000);

    // // Print the resulting vector y
    // PetscPrintf(PETSC_COMM_SELF, "Resulting vector y:\n");
    // VecView(y, PETSC_VIEWER_STDOUT_SELF);

    // // Clean up
    // MatDestroy(&A);
    // MatDestroy(&B);
    // MatDestroy(&K);
    // VecDestroy(&x);
    // VecDestroy(&y);

    // PetscFinalize();
    // return 0;
}
