#include "particle_list.hpp"
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"
#include "../grid/grid.hpp"

Particle_Field_List::Particle_Field_List(long np){
    np_=np;

	parts_.reserve((long)1.5*np); // Have at least 1.5x the number of particles for 
	                      // slosh room. Excessive maybe?
}

Particle_Field_List::~Particle_Field_List(){
}

void Particle_Field_List::Load(int restart){
    //dummy code inserted by Yuan for testing main.cpp
    if(restart==0){// initial run
       for(long ip=0;ip<np_;ip++){
           Particle* p = new_particle();
           p->x1=-1.0;
           p->x2=1.0;
           p->x3=ip*1.0;
    
           p->v1=0.0;
           p->v2=4.0;
           p->v3=1.0;
           parts_.push_back(p);
       }
    } 
}

void Particle_Field_List::Push(double dt){

    for(long ip=0;ip<np_;ip++){
        pusher_->Step(parts_[ip],parts_[ip]->field,dt);
    }

}

void Particle_Field_List::Pass(){
}

long Particle_Field_List::nParticles(){
	return np_;
}

void Particle_Field_List::InterpolateEB(Grid* grid){
  Interpolator *interpolator = new Interpolator();

  long iCell = 0; //cell # tracker
  long pCell = 0; //particle cell #
  double cellvars[21];//Vector describing position of and all field elements of a cell
                      //["least" corner vertex, E-field on edges, B-field on surfaces]
                      //21 elements ordered as: [x1,x2,x3,E1,E2,E3,B1,B2,B3]
                      //              of sizes:   1, 1, 1, 4, 4, 4, 2, 2, 2
  double pos[3]; //Vector of position of particle.
  double lcell[3]; //Vector of lengths of cell.

  //Get lengths of grid cells.
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  for (long i=0; i<np_; i++) {
    //Get position of particle.
    pos[0] = parts_[i]->x1;
    pos[1] = parts_[i]->x2;
    pos[2] = parts_[i]->x3;

    //Update cell field variables.
    pCell = grid->getCellID(pos[0],pos[1],pos[2]);
    if (pCell >= 0) {
      if (pCell != iCell) {
	iCell = pCell;
	grid->getFieldInterpolatorVec(iCell, cellvars);
      }

      //Interpolate fields at particle.
      interpolator->interpolate_fields(pos, lcell, cellvars, parts_[i]->field);
    }
  }
}

void Particle_Field_List::SortParticles(Particle_Compare comp){
	std::sort(parts_.begin(),parts_.end(),comp);
}

void Particle_Field_List::depositRhoJ(Grid *grid, double dt){
  Depositor *depositor = new Depositor();

  long cellID[2] = {-1, -1}; //cell id tracker, 2 slots: 'exited' cell and 'entered' cell
  long tempCellID[2] = {}; //particle's 'exited' and 'entered' cell ids.
  double cvec[3] = {};
  double cellverts[2][3] = {};//Vector describing two cells' respective "least" vertices.
  double pos[2][3] = {}; //Array of previous and current positions of particle.
  double dpos[3] = {}; //Vector of change in position of particle.
  double lcell[3] = {}; //Vector of lengths of cells.
  double RhoJObj[2][12] = {}; //Array describing two cells' Jx (4), Jy (4), Jz (4).
  double pcharge;

  //Get lengths of grid cells.
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  //Zero the grid's currents and charge densities.
  grid->zeroRhoJ();

  //Cycle through particles, depositing RhoJ for each one.
  for (long i=0; i<np_; i++) {
    //Get last and current positions of particle.
    pos[0][0] = parts_[i]->xo1;
    pos[0][1] = parts_[i]->xo2;
    pos[0][2] = parts_[i]->xo3;
    pos[1][0] = parts_[i]->x1;
    pos[1][1] = parts_[i]->x2;
    pos[1][2] = parts_[i]->x3;
    dpos[0] = parts_[i]->dx1;
    dpos[1] = parts_[i]->dx2;
    dpos[2] = parts_[i]->dx3;

    //Get charge of particle.
    pcharge = parts_[i]->q;

    //Update particle's last and current cell id's
    for (int j=0; j<2; j++) tempCellID[j] = grid->getCellID(pos[j][0],pos[j][1],pos[j][2]);

    //Only if particle ENDS in a 'real' (non-ghost) cell, do we 'deposit' it.
    if (tempCellID[1] >= 0) {
      //If list of particles continue to start and end in same cells, add current to existing RhoJObj object.
      //Otherwise, deposit (add) RhoJ to the grid, and re-point (and zero out) the existing RhoJObj object.
      for (int j=0; j<2; j++) {
	if (tempCellID[j] != cellID[j]) {
	  //If cellID has already been assinged...
	  if (cellID[j] != -1) {
	    //Deposit to grid
	    grid->addRhoJ(cellID[j],&RhoJObj[j][0]);
	    //Zero RhoJObj[j]
	    for (int k=0; k<12; k++) {
	      RhoJObj[j][k] = 0;
	    }
	  }
	  //Get new cell vertex data.
	  cellID[j] = tempCellID[j];
	  grid->getCellVertex(cellID[j], cvec);
	  for (int k=0; k<3; k++) {
	    cellverts[j][k] = cvec[k];
	  }
	}
      }

      //Generate currents at cell edges.
      depositor->deposit_particle_RhoJ(cellID, pos, dpos, lcell, cellverts, dt, pcharge, RhoJObj);
    }
  }

  //Add remaining currents to the grid.
  for (int j=0; j<2; j++) {
    grid->addRhoJ(cellID[j],&RhoJObj[j][0]);
  }
}

double Particle_Field_List::maxVelocity(void){
    return 0.1;
}
