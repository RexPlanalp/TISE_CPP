#include "tise.h"
#include <petscsys.h>
#include <slepceps.h>
#include <vector>
#include "basis.h"
#include <mpi.h>
#include <complex>

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    Mat K;
    Mat V;
    Mat S;

    ierr = SlepcInitialize(&argc, &argv, NULL,NULL); CHKERRQ(ierr);

    int n_basis = 1500;
    int degree = 6;
    double rmax = 1500;
    int lmax = 75;
    int nmax = 10;
    double R0 = 1245;
    double eta=  0;
    double tise_tolerance = 1E-15;
    PetscInt tise_mat_iter = 3000;

    // int n_basis = 3;
    // int degree = 2;
    // double rmax = 3;
    // int lmax = 1;
    // int nmax = 10;
    // double R0 = 3;
    // double eta=  0;
    // double tise_tolerance = 1E-15;
    // PetscInt tise_mat_iter = 3000;




    bool SAVE_BSPLINES = false;
    bool EMBED = false;


    std::vector<PetscScalar> knots = basis::linear_knots(n_basis, degree, rmax);
    if (SAVE_BSPLINES)
    {
        int Nx = rmax * 10;
        int rex_err = basis::save_bsplinee_basis(n_basis,degree,knots,Nx,rmax,R0,eta);
    }
  
    PetscPrintf(PETSC_COMM_WORLD, "Solving TISE\n");
    double start_time = MPI_Wtime();
    
    tise::solve_tise(n_basis,degree,knots,lmax,nmax,R0,eta,tise_tolerance,tise_mat_iter,EMBED);

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

    ierr = SlepcFinalize(); CHKERRQ(ierr);
    return 0;
}
