#include <stdlib.h>
#include <math.h>
#include "poisson.hpp"

Poisson_Solver::Poisson_Solver(int *nxyz, int nGhosts, double *xyz0, double *Lxyz)
  : Grid(nxyz, nGhosts, xyz0, Lxyz) {

  int ifield = -1;
  rho_=newField_(++ifield);
  phi1_=newField_(++ifield);
  phi2_=newField_(++ifield);
  Ax1_=newField_(++ifield);
  Ay1_=newField_(++ifield);
  Az1_=newField_(++ifield);
  Ax2_=newField_(++ifield);
  Ay2_=newField_(++ifield);
  Az2_=newField_(++ifield);
}

void Poisson_Solver::run_poisson_solver(double*** u0, double*** u1,double*** R,double convergenceTol) {
  //u0 is the first guess at a solution to Poisson's equation.
  //u1 is the work array of equal size
  //R is the 'source' array for Poisson's eqn

  double lcell[3] = {}; //Vector of lengths of cells.                                                                             

  double ncells[3] = {};
  int err = getNumCells3D(ncells);
  long nx = ncells[0];
  long ny = ncells[1];
  long nz = ncells[2];

  for (int i=0; i<3; i++) lcell[i] = getStepSize(i);

  //Define constants used to iterate Poisson's equation
  double celldist2 = pow(lcell[0], 2.0) + pow(lcell[1], 2.0) + pow(lcell[2], 2.0);
  double ax = pow(lcell[1], 2.0) * pow(lcell[2], 2.0) / (2.0 * celldist2);
  double ay = pow(lcell[0], 2.0) * pow(lcell[2], 2.0) / (2.0 * celldist2);
  double az = pow(lcell[0], 2.0) * pow(lcell[1], 2.0) / (2.0 * celldist2);
  double af = pow(lcell[0], 2.0) * pow(lcell[1], 2.0) * pow(lcell[2], 2.0) / (2.0 * celldist2);

  //Understand how to get and use the u0, u1 and Rho arrays used below...

  bool jacobi_method_converged = false;
  double maxDiff = 0.0;
  double absDiff = 0.0;

  int iternum = -1;
  do {
    iternum++;
    //iterate over entire grid. Note boundary conditions must be supplied!
    for ( int i=1; i<nx; i++ ) {
      for ( int j=1; j<ny; j++ ) {
	for ( int k=1; k<nz; k++ ) {
	  if ( iternum % 2 == 0 ) {
	    u1[i][j][k] = ax*(u0[i-1][j][k]+u0[i+1][j][k]) + ay*(u0[i][j-1][k]+u0[i][j+1][k]) + 
	      az*(u0[i][j][k-1]+u0[i][j][k+1]) - af*R[i][j][k];
	  } else {
	    u0[i][j][k] = ax*(u1[i-1][j][k]+u1[i+1][j][k]) + ay*(u1[i][j-1][k]+u1[i][j+1][k]) +
              az*(u1[i][j][k-1]+u1[i][j][k+1]) - af*R[i][j][k];
	  }

	  absDiff = abs(u1[i][j][k] - u0[i][j][k]);
	  if (absDiff > maxDiff) maxDiff = absDiff;
	}
      }
    }

    //Pass and receive MPI for boundary!  Pass the boundary first and Receive the boundary last!
    //Determine GLOBAL convergence of jacobi method across all MPI domains!!!

    if (maxDiff < convergenceTol) jacobi_method_converged = true;
  }while( !jacobi_method_converged );

  //If last calculated pass updated the work array (u1), copy it to the solution array (u0).
  if ( iternum % 2 == 0 ) {
    for ( int i=1; i<nx; i++ ) {
      for ( int j=1; j<ny; j++ ) {
        for ( int k=1; k<nz; k++ ) {
	  u0[i][j][k] = u1[i][j][k];
	}
      }
    }
  }

}
