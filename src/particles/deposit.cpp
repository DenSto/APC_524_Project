#include <stdlib.h>
#include "deposit.hpp"

Depositor::Depositor(){
}

Depositor::~Depositor(){
}

void Depositor::deposit_particle_RhoJ(Particle *part, double* lcell, double* cellverts, double* RhoJObj) {
  //Function returns an updated RhoJObj to include the current due to particle motion on the cell edges.

  double curIn1D = 0.0; //1D current
  double cur2AreaRatio = 0.0; //Current density to area ratio.
  double cellArea = 0.0; //Cell face area.
  double cellVolume = lcell[0]*lcell[1]*lcell[2]; //Cell volume.
  double relpos[3] = {}; //Relative position of particle to 'least' vertex of cell.

  for (int i=0; i<3; i++) relpos[i] = part->x[i] - cellverts[i];

  for (int i=0; i<3; i++) {
    int j, k;
    j = i+1 % 3; /*Modular arithmetic cycles over Cartesian directions*/
    k = i+2 % 3;

    //Calculate current due to instantaneous velocity and apply to entered (current) cell.
    cellArea = lcell[j]*lcell[k];
    curIn1D = part->q * part->v[i];
    cur2AreaRatio = (curIn1D / cellVolume) / cellArea;
    RhoJObj[4*i] += cur2AreaRatio * (lcell[j]-relpos[j]) * (lcell[k]-relpos[k]);
    RhoJObj[4*i+1] += cur2AreaRatio * relpos[j] * (lcell[k]-relpos[k]);
    RhoJObj[4*i+2] += cur2AreaRatio * relpos[j] * relpos[k];
    RhoJObj[4*i+3] += cur2AreaRatio * (lcell[j]-relpos[j]) * relpos[k];
  }
}
