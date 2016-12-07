#include "particle_list.hpp"
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"
#include "../grid/grid.hpp"

Particle_List::Particle_List(long np){
    np_=np;

	parts_.reserve((long)1.5*np); // Have at least 1.5x the number of particles for 
	                      // slosh room. Excessive maybe?
}

Particle_List::~Particle_List(){
}

void Particle_List::Load(int restart){
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

void Particle_List::Push(double dt){

    for(long ip=0;ip<np_;ip++){
        pusher_->Step(parts_[ip],parts_[ip]->field,dt);
    }

}

void Particle_List::Pass(){
}

long Particle_List::nParticles(){
	return np_;
}

void Particle_List::InterpolateEB(Grid* grid){
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
    if (pCell != iCell) {
      iCell = pCell;
      grid->getFieldInterpolatorVec(iCell, cellvars);
    }

    //Interpolate fields at particle.
    interpolator->interpolate_fields(pos, lcell, cellvars, parts_[i]->field);
  }
}

void Particle_List::SortParticles(Particle_Compare comp){
	std::sort(parts_.begin(),parts_.end(),comp);
}

void Particle_List::depositCurrent(Grid *grids){
}

void Particle_List::depositCharge(Grid *grids){
}

double Particle_List::maxVelocity(void){
    return 0.1;
}
