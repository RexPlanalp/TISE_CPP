#include "tise.h"
#include <petscsys.h>
#include <slepceps.h>
#include <vector>
#include "basis.h"
#include <mpi.h>
#include <complex>

PetscInt main(PetscInt argc, char **argv) {
    
    PetscErrorCode ierr;
    ierr = SlepcInitialize(&argc, &argv, NULL,NULL); CHKERRQ(ierr);

    PetscInt n_basis = 3000;
    PetscInt degree = 6;
    double rmax = 3000;
    PetscInt lmax = 75;
    PetscInt nmax = 10;
    double R0 = 3000;
    double eta=  0;
    double tise_tolerance = 1E-15;
    PetscInt  tise_mat_iter = 3000;
    bool SAVE_BSPLINES = false;

    std::vector<PetscScalar> knots = basis::linear_knots(n_basis, degree, rmax);

    if (SAVE_BSPLINES)
    {
        PetscInt Nx = rmax * 10;
        ierr = basis::save_bsplinee_basis(n_basis,degree,knots,Nx,rmax,R0,eta); CHKERRQ(ierr);
    }
  
    PetscPrintf(PETSC_COMM_WORLD, "Solving TISE\n");
    double start_time = MPI_Wtime();
    tise::solve_tise(n_basis,degree,knots,lmax,nmax,R0,eta,tise_tolerance,tise_mat_iter);
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

    ierr = SlepcFinalize(); CHKERRQ(ierr);
    return 0;
}
