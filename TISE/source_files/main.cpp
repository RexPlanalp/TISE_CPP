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

    // int n_basis = 1250;
    // int degree = 6;
    // double rmax = 1250;
    // int lmax = 50;
    // int nmax = 10;
    // double R0 = 1245;
    // double eta=  0;

    //For testing
    // int n_basis = 30;
    // int degree = 6;
    // double rmax = 10;
    // int lmax = 5;
    // int nmax = 5;
    // double R0 = 5;
    // double eta=  0;


    int n_basis = 4;
    int degree = 2;
    double rmax = 10;
    int lmax = 5;
    int nmax = 5;
    double R0 = 5;
    double eta=  0.25;

    std::vector<std::complex<double>> knots = basis::linear_knots(n_basis, degree, rmax);


    // std::complex<double> result = basis::overlap_matrix_element(3, 3, degree, knots, R0, eta);
    // std::cout << result << std::endl;



    // double x_test = 6;
    // std::complex<double> x_modified = basis::R_x(x_test,R0,eta);
    // std::vector<std::complex<double>> knots_modified = basis::R_k(knots,R0,eta);
    // std::complex<double> result = basis::B(3,degree,knots_modified,x_modified);

    // std::cout << result << std::endl;






    // For testing
    // int Nx = 1000;
    // int rex_err = basis::save_bsplinee_basis(n_basis,degree,knots,Nx,rmax,R0,eta);



    PetscPrintf(PETSC_COMM_WORLD, "Solving TISE\n");


    double start_time = MPI_Wtime();
    
    tise::solve_tise(n_basis,degree,knots,lmax,nmax,R0,eta);


    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    PetscPrintf(PETSC_COMM_WORLD, "Time taken to construct matrix: %.3f seconds\n", elapsed_time);

    
    ierr = SlepcFinalize(); CHKERRQ(ierr);

    return 0;
}
