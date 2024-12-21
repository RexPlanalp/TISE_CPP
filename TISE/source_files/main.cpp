#include "tise.h"
#include <petscsys.h>
#include <vector>
#include "basis.h"
#include <mpi.h>

int main(int argc, char **argv) {
    

    PetscErrorCode ierr;
    Mat K;
    Mat V;
    Mat S;

    int n_basis = 3000;
    int degree = 6;
    double rmax = 3000;

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
    std::vector<double> knots = basis::linear_knots(n_basis, degree, rmax);

    double start_time = MPI_Wtime();
    ierr = tise::construct_kinetic_matrix(&K, n_basis, degree, knots); CHKERRQ(ierr);
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

    start_time = MPI_Wtime();
    ierr = tise::construct_kinetic_matrix(&V, n_basis, degree, knots); CHKERRQ(ierr);
    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

    start_time = MPI_Wtime();
    ierr = tise::construct_overlap_matrix(&S, n_basis, degree, knots); CHKERRQ(ierr);
    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);


    // Clean up and finalize
    ierr = MatDestroy(&K); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
