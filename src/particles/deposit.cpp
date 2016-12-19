#include <stdlib.h>
#include "deposit.hpp"

Depositor::Depositor(){
}

Depositor::~Depositor(){
}

void Depositor::deposit_particle_RhoJ(Particle *part, double* lcell, double* cellverts, double* RhoJObj) {
  //Function returns an updated RhoJObj to include the current due to particle motion on the cell edges.

  double curIn1D = 0.0;
  double cur2AreaRatio = 0.0;
  double cellArea = 0.0;

  for (int i=0; i<3; i++) {
    int j, k;
    j = i+1 % 3; /*Modular arithmetic cycles over Cartesian directions*/
    k = i+2 % 3;

    //Calculate current due to instantaneous velocity and apply to entered (current) cell.
    cellArea = lcell[j]*lcell[k];
    curIn1D = part->q * part->v[i];
    cur2AreaRatio = curIn1D / cellArea;
    RhoJObj[4*i] += cur2AreaRatio * (lcell[j]-part->x[j]) * (lcell[k]-part->x[k]);
    RhoJObj[4*i+1] += cur2AreaRatio * part->x[j] * (lcell[k]-part->x[k]);
    RhoJObj[4*i+2] += cur2AreaRatio * part->x[j] * part->x[k];
    RhoJObj[4*i+3] += cur2AreaRatio * (lcell[j]-part->x[j]) * part->x[k];
  }
}
