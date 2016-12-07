#include "particle_handler.hpp"
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"
#include "../grid/grid.hpp"

Particle_Handler::Particle_Handler(long np){
    np_=np;

	parts_.reserve((long)1.5*np); // Have at least 1.5x the number of particles for 
	                      // slosh room. Excessive maybe?
}

Particle_Handler::~Particle_Handler(){
}

void Particle_Handler::Load(int restart){
    //dummy code inserted by Yuan for testing main.cpp
    if(restart==0){// initial run
       for(long ip=0;ip<np_;ip++){
           Particle p = new_particle();
           p.x[0]=-1.0;
           p.x[1]=1.0;
           p.x[2]=ip*1.0;
    
           p.v[0]=0.0;
           p.v[1]=4.0;
           p.v[2]=1.0;
           parts_.push_back(p);
       }
    } 
}

void Particle_Handler::Push(double dt){

    for(long ip=0;ip<np_;ip++){
        pusher_->Step(&(parts_[ip]),&(parts_[ip].field),dt);
    }

}

void Particle_Handler::Pass(){
}

long Particle_Handler::nParticles(){
	return np_;
}

void Particle_Handler::InterpolateEB(Grid* grid){
  Interpolator *interpolator = new Interpolator();

  long iCell = 0; //cell # tracker
  long pCell = 0; //particle cell #
  double cellvars[21];//Vector describing position of and all field elements of a cell
                      //["least" corner vertex, E-field on edges, B-field on surfaces]
                      //21 elements ordered as: [x[0],x[1],x[2],E1,E2,E3,B1,B2,B3]
                      //              of sizes:   1, 1, 1, 4, 4, 4, 2, 2, 2
  double pos[3]; //Vector of position of particle.
  double lcell[3]; //Vector of lengths of cell.

  //Get lengths of grid cells.
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  for (long i=0; i<np_; i++) {
    //Get position of particle.
    pos[0] = parts_[i].x[0];
    pos[1] = parts_[i].x[1];
    pos[2] = parts_[i].x[2];

    //Update cell field variables.
    pCell = grid->getCellID(pos[0],pos[1],pos[2]);
    if (pCell != iCell) {
      iCell = pCell;
      grid->getFieldInterpolatorVec(iCell, cellvars);
    }

    //Interpolate fields at particle.
    interpolator->interpolate_fields(pos, lcell, cellvars, &(parts_[i].field));
  }
}

void Particle_Handler::SortParticles(Particle_Compare comp){
	std::sort(parts_.begin(),parts_.end(),comp);
}

void Particle_Handler::depositCurrent(Grid *grids){
}

void Particle_Handler::depositCharge(Grid *grids){
}

double Particle_Handler::maxVelocity(void){
    return 0.1;
}
