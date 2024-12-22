#include "tise.h"
#include <petscsys.h>
#include <slepceps.h>
#include <vector>
#include "basis.h"
#include <mpi.h>

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    Mat K;
    Mat V;
    Mat S;

    int n_basis = 4;
    int degree = 2;
    double rmax = 4;
    int lmax = 2;
    int nmax = 2;

    




    ierr = SlepcInitialize(&argc, &argv, NULL,NULL); CHKERRQ(ierr);
    std::vector<double> knots = basis::linear_knots(n_basis, degree, rmax);

    PetscPrintf(PETSC_COMM_WORLD, "Generating Bsplines\n");


    double start_time = MPI_Wtime();
    
    tise::solve_tise(n_basis,degree,knots,lmax,nmax);


    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

    
    ierr = SlepcFinalize(); CHKERRQ(ierr);

    return 0;
}
