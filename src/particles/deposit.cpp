#include <stdlib.h>
#include <stdio.h>
#include "deposit.hpp"
#include "../globals.hpp"

Depositor::Depositor(){
}

Depositor::~Depositor(){
}

void Depositor::deposit_particle_J(Particle *part, double* lcell, double* cellverts, double* JObj) {
  //Function returns an updated JObj to include the current due to particle motion on the cell edges.

  double pcharge = 0.0; //particle charge
  double curIn1D = 0.0; //1D current
  double cur2AreaRatio = 0.0; //Current density to area ratio.
  double cellArea = 0.0; //Cell face area.
  double cellVolume = lcell[0]*lcell[1]*lcell[2]; //unit cell volume.
  double relpos[3] = {}; //Relative position of particle to 'least' vertex of cell.

  for (int i=0; i<3; i++) relpos[i] = part->x[i] - cellverts[i];

  pcharge = part->q;
  for (int i=0; i<3; i++) {
    int j, k;
    j = (i+1) % 3; /*Modular arithmetic cycles over Cartesian directions*/
    k = (i+2) % 3;

    //Calculate current due to instantaneous velocity and apply to entered (current) cell.
    cellArea = lcell[j]*lcell[k];
    curIn1D = UNIT_RHOJ * pcharge * part->v[i]; // multiply unit of current
    cur2AreaRatio = (curIn1D / cellVolume) / cellArea;
    JObj[4*i] += cur2AreaRatio * (lcell[j]-relpos[j]) * (lcell[k]-relpos[k]);
    JObj[4*i+1] += cur2AreaRatio * relpos[j] * (lcell[k]-relpos[k]);
    JObj[4*i+2] += cur2AreaRatio * relpos[j] * relpos[k];
    JObj[4*i+3] += cur2AreaRatio * (lcell[j]-relpos[j]) * relpos[k];
  }
}

void Depositor::deposit_particle_Rho(Particle *part, double* lcell, double* cellverts, double* RhoObj) {
  //Function returns an updated RhoObj to include the charge due to particles on the cell edges.

  double pcharge = 0.0; //particle charge
  double charge2VolRatio = 0.0; //Current density to volume ratio.
  double cellVolume = lcell[0]*lcell[1]*lcell[2]; //Cell volume.
  double relpos[3] = {}; //Relative position of particle to 'least' vertex of cell.

  for (int i=0; i<3; i++) relpos[i] = part->x[i] - cellverts[i];

  //Calculate charge and apply to entered (current) cell.
  pcharge = part->q;
  charge2VolRatio = UNIT_RHOJ*(pcharge / cellVolume) / cellVolume; //multiply unit of Rho
  RhoObj[0] += charge2VolRatio * (lcell[0]-relpos[0]) * (lcell[1]-relpos[1]) * (lcell[2]-relpos[2]); //[0 0 0]
  RhoObj[1] += charge2VolRatio * relpos[0] * (lcell[1]-relpos[1]) * (lcell[2]-relpos[2]); //[1 0 0]
  RhoObj[2] += charge2VolRatio * relpos[0] * relpos[1] * (lcell[2]-relpos[2]); //[1 1 0]
  RhoObj[3] += charge2VolRatio * (lcell[0]-relpos[0]) * relpos[1] * (lcell[2]-relpos[2]); //[0 1 0]
  RhoObj[4] += charge2VolRatio * (lcell[0]-relpos[0]) * relpos[1] * relpos[2]; //[0 1 1]
  RhoObj[5] += charge2VolRatio * relpos[0] * relpos[1] * relpos[2]; //[1 1 1]
  RhoObj[6] += charge2VolRatio * relpos[0] * (lcell[1]-relpos[1]) * relpos[2]; //[1 0 1]
  RhoObj[7] += charge2VolRatio * (lcell[0]-relpos[0]) * (lcell[1]-relpos[1]) * relpos[2]; //[0 0 1]
}
