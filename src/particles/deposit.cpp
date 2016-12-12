#include <stdlib.h>
#include "deposit.hpp"
#include "assert.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

Depositor::Depositor(){
}

Depositor::~Depositor(){
}

//Comment(YShi): J deposition is better done with velocity directly. Because at 
//initialization step, dpos is undefined. The interface can also be made simpler
//by passing in struct Particle, which contains pos, dpos,q 
//If the above is accepted, then some contents of this soutine needs to be changed
void Depositor::deposit_particle_RhoJ(long* cellID,Particle *part, double* lcell, double cellverts[][3], double RhoJObj[][12]) {
//void Depositor::deposit_particle_RhoJ(long* cellID, double pos[][3], double* dpos, double* lcell, double cellverts[][3], double dt, double q, double RhoJObj[][12]) {
  //Function returns an updated RhoJObj to include the current due to particle motion on the cell edges.

//  double dx = 0.0;
//  double fullArea = 0.0;
  double curIn1D = 0.0;
  double cur2AreaRatio = 0.0;
  double avgpos[3] = {};
  int celltot = 0;
  double jFactor = 0.0;
  double moveDir = 0.0;
  double cellArea = 0.0;

  //If particle starts and ends in the same cell.
  if (cellID[0]==cellID[1]) {
    celltot = 1;
    jFactor = 1.0;
  } else {
    celltot = 2;
    jFactor = 0.5;
  }

  assert(celltot > 0 && celltot <= 2); //Other cases are not yet implemented!
  for (int kcell=0; kcell<celltot; kcell++) {
    for (int i=0; i<3; i++) {
      int j, k;
      j = i+1 % 3; /*Modular arithmetic cycles over Cartesian directions*/
      k = i+2 % 3;

      //Get direction of particle movement w.r.t. point of calculation. (Propagates backward from endpoint, e.g.)
      if (kcell == 0) {
	moveDir = 1.0;
      } else  if (kcell == 1) {
	moveDir = -1.0;
      }

      //Get average position of particle w.r.t. "least" cell vertex.
/*      for (int j=0; j<3; j++) {
	avgpos[j] = ( pos[kcell][j] + (pos[kcell][j] + moveDir*dpos[j]) )/2.0 - cellverts[kcell][j];
	//Keep its propagation within the cell, however.
	avgpos[j] = MIN(MAX(avgpos[j], 0.0), lcell[j]);
      }
*/
      //Calculate current DUE TO MOTION and apply to exited and entered cells (split in half if not the same).
      cellArea = lcell[j]*lcell[k];
//      curIn1D = jFactor*q*dpos[i]/dt;
      cur2AreaRatio = curIn1D / cellArea;
      RhoJObj[kcell][4*i] += cur2AreaRatio * (lcell[j]-avgpos[j]) * (lcell[k]-avgpos[k]);
      RhoJObj[kcell][4*i+1] += cur2AreaRatio * avgpos[j] * (lcell[k]-avgpos[k]);
      RhoJObj[kcell][4*i+2] += cur2AreaRatio * avgpos[j] * avgpos[k];
      RhoJObj[kcell][4*i+3] += cur2AreaRatio * (lcell[j]-avgpos[j]) * avgpos[k];
    }
  }
}
