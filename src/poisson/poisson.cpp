#include <stdlib.h>
#include <math.h>
#include "poisson.hpp"

Poisson_Solver::Poisson_Solver(){
}

Poisson_Solver::~Poisson_Solver(){
}

void Poisson_Solver::run_poisson_solver(Grid* grid, Domain* domain) {

  double lcell[3] = {}; //Vector of lengths of cells.                                                                             

  double ncells[3] = {};
  int err = grid->getNumCells3D(ncells);
  long nx = ncells[0];
  long ny = ncells[1];
  long nz = ncells[2];

  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  //Define constants used to iterate Poisson's equation
  double celldist2 = pow(lcell[0], 2.0) + pow(lcell[1], 2.0) + pow(lcell[2], 2.0);
  double ax = pow(lcell[1], 2.0) * pow(lcell[2], 2.0) / (2.0 * celldist2);
  double ay = pow(lcell[0], 2.0) * pow(lcell[2], 2.0) / (2.0 * celldist2);
  double az = pow(lcell[0], 2.0) * pow(lcell[1], 2.0) / (2.0 * celldist2);
  double af = pow(lcell[0], 2.0) * pow(lcell[1], 2.0) * pow(lcell[2], 2.0) / (2.0 * celldist2);

  //Understand how to get and use the u0, u1 and Rho arrays used below...


  bool jacobi_method_converged = false;
  int iternum = 0;
  do {
    //iterate over entire grid. Note boundary conditions must be supplied!
    for ( int i=1; i<nx; i++ ) {
      for ( int j=1; j<ny; j++ ) {
	for ( int k=1; k<nz; k++ ) {
	  if ( iternum % 2 == 0 ) {

	    //	    u1[i][j][k] = ax*(u0[i-1][j][k]+u0[i+1][j][k]) + ay*(u0[i][j-1][k]+u0[i][j+1][k]) + 
	    //	      az*(u0[i][j][k-1]+u0[i][j][k+1]) - af*Rho[i][j][k];

	  } else {
	    //Once the above is finalized...implement the same here assiging u1 to u0...

	    //u0[i][j][k] = ax*(u1[i-1][j][k]+u1[i+1][j][k]) + ay*(u1[i][j-1][k]+u1[i][j+1][k]) +
            //  az*(u1[i][j][k-1]+u1[i][j][k+1]) - af*Rho[i][j][k];
	  }
	}
      }
    }

    //Pass and receive MPI for boundary!  Pass the boundary first and Receive the boundary last!
    //Determine GLOBAL convergence of jacobi method across all MPI domains!!!
  }while( !jacobi_method_converged );

}
