#include "../globals.hpp"
#include "particle_handler.hpp"
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"
#include "../grid/grid.hpp"
#include "../utils/RNG.hpp"
#define _USE_MATH_DEFINES
#if USE_MPI
#include "mpi.h"
#endif

Particle_Handler::Particle_Handler(){
	np_=0;
	//parts_.reserve((long)1.5*np); // Have at least 1.5x the number of particles for 
                     // slosh room. Excessive maybe?
}

Particle_Handler::~Particle_Handler(){
}

void Particle_Handler::Load(Input_Info_t *input_info, Domain* domain){

    int restart = input_info->restart;
    double* L = domain->getLxyz();
    double* x0 = domain->getxyz0();

    int nspec      = input_info->nspecies;
    int npart      = input_info->np;
    double *mass   = input_info->mass_ratio;
    double *charge = input_info->charge_ratio;
    double *dens   = input_info->dens_frac; 

    if(restart==0){// initial run
	Random_Number_Generator *rng = new Random_Number_Generator(-1);
        int ispec = 0; // temporaty counter
        double cden = dens[0]; // cummulative density fraction
        double vth;
        for(long ip=0; ip < npart;ip++){
		Particle p = new_particle();
                if(ip >= cden*npart){
			ispec += 1;
                        cden  += dens[ispec];
		}
		p.q = charge[ispec];
		p.m = mass[ispec];

		//double vth= sqrt(8 * input_info->temp / (p.m * M_PI));
		vth=UNIT_VTH*sqrt(input_info->temp[ispec]/p.m);

		p.x[0]=rng->getUniform()*L[0]+x0[0];
		p.x[1]=rng->getUniform()*L[1]+x0[1];
		p.x[2]=rng->getUniform()*L[2]+x0[2];

		p.v[0]=rng->getGaussian(0.0,vth);
		p.v[1]=rng->getGaussian(0.0,vth);
		p.v[2]=rng->getGaussian(0.0,vth);

		parts_.push_back(p);
		np_++;
    } }
    else{//read restart file
        //dummy code inserted by YShi for testing
        //insert a single particle at the center of the cell
	Particle p = new_particle();
	p.q = 0;
	p.m = 1;
	p.x[0]=L[0]/2+x0[0];
	p.x[1]=L[1]/2+x0[1];
	p.x[2]=L[2]/2+x0[2];
//	p.v[0]=0.1;
	p.v[1]=0.1;
//	p.v[2]=0.1;
	parts_.push_back(p);
	np_++;
    }
}

void Particle_Handler::Push(double dt){

    for(long ip=0;ip<np_;ip++){
        pusher_->Step(&(parts_[ip]),&(parts_[ip].field),dt);
    }

}

long Particle_Handler::nParticles(){
	return np_;
}

void Particle_Handler::incrementNParticles(int inc){
	np_+=inc;
}

void Particle_Handler::InterpolateEB(Grid* grid){
  fprintf(stderr,"rank=%d,start InterpolateEB\n",rank_MPI);
  Interpolator *interpolator = new Interpolator();
  fprintf(stderr,"rank=%d,new Interpolator\n",rank_MPI);

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
    //fprintf(stderr,"rank=%d,pCell=%ld\n",rank_MPI,pCell);
    if (pCell >= 0) {
      if (pCell != iCell) {
	iCell = pCell;
	grid->getFieldInterpolatorVec(iCell, cellvars);
      }

      //Interpolate fields at particle.
      fprintf(stderr,"rank=%d,call Interpolator\n",rank_MPI);
      interpolator->interpolate_fields(pos, lcell, cellvars, &(parts_[i].field));
    }
  }
  fprintf(stderr,"rank=%d,Finish InterpolateEB\n",rank_MPI);
}

void Particle_Handler::SortParticles(Particle_Compare comp){
	std::sort(parts_.begin(),parts_.end(),comp);
}

void Particle_Handler::depositRhoJ(Grid *grid, bool depositRho){
  Depositor *depositor = new Depositor();

  long cellID = -1; //cell id tracker.
  long tempCellID = 0; //particle's 'entered' cell id.
  double cellverts[3] = {};//Vector describing a cell's "least" vertex.
  double lcell[3] = {}; //Vector of lengths of cells.
  double JObj[12] = {}; //Array describing a cell's Jx (4), Jy (4), Jz (4).
  double RhoObj[8] = {}; //Array describing a cell's 8 Rho vertices.

  //Get lengths of grid cells.
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  //Zero the grid's currents and charge densities.
  grid->zeroJ();
  if (depositRho) grid->zeroRho();

  //Cycle through particles, depositing RhoJ for each one.
  for (long i=0; i<np_; i++) {
    //Get particle's 'entered' (post-push) cell id
    tempCellID = grid->getCellID(parts_[i].x[0],parts_[i].x[1],parts_[i].x[2]);

    //Only if particle ENDS in a 'real' (non-ghost) cell, do we 'deposit' it.
    if (tempCellID >= 0) {
      //If list of particles continues to be in same cell, add current to existing JObj object.
      //Otherwise, deposit (add) RhoJ to the grid, and re-point (and zero out) the existing JObj object.
      if (tempCellID != cellID) {
	//If cellID has already been assigned...
	if (cellID != -1) {
	  //Deposit to grid
	  grid->addJ(cellID,&JObj[0]);
	  if (depositRho) grid->addRho(cellID,&RhoObj[0]);
	  //Zero JObj
	  for (int k=0; k<12; k++) JObj[k] = 0;
	  if (depositRho) {
	    for (int k=0; k<8; k++) RhoObj[k] = 0;
	  }
	}
	//Get new cell vertex data.
	cellID = tempCellID;
	grid->getCellVertex(cellID, cellverts);
      }

      //Generate currents at cell edges.
      depositor->deposit_particle_J(&(parts_[i]), lcell, cellverts, JObj);
      if (depositRho) depositor->deposit_particle_Rho(&(parts_[i]), lcell, cellverts, RhoObj);
    }
  }

  //Add remaining current to the grid.
  grid->addJ(cellID,&JObj[0]);
  if (depositRho) grid->addRho(cellID,&RhoObj[0]);
}

double Particle_Handler::computeCFLTimestep(Domain* domain){
/*	double maxV[3], mindt[3];

	double *dx = domain->getdx();
	for(int i = 0; i < 3; i++) maxV[i]=0.0;
	for(std::vector<Particle>::iterator iter = parts_.begin(); iter != parts_.end(); iter++){
		for(int i = 0; i < 3; i++){
			if(fabs(iter->v[i]) > maxV[i])		
				maxV[i] = fabs(iter->v[i]);
		}			
	}
	for(int i = 0; i < 3; i++){
		assert(maxV[i] > 0);
		mindt[i] = dx[i]/maxV[i];
	}
	delete[] dx;
	
#if USE_MPI
//	double mindtall[3];
//	int ierr = MPI_Allreduce(mindt,mindtall,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	//return *std::min_element(mindtall,mindtall+3);
	return *std::min_element(mindt,mindt+3);
#else
	return *std::min_element(mindt,mindt+3);
#endif
*/
	return 0.01;
}


void Particle_Handler::clearGhosts(){
	for(std::vector<Particle>::iterator iter = parts_.begin(); iter != parts_.end();){
		if(iter->isGhost){
			std::swap(*iter, parts_.back());
			parts_.pop_back();
		} else {
			iter++;
		}
	}
	assert(parts_.size() == np_);
//        std::cerr<<"parts_.size="<<parts_.size()<<", np_="<<np_<<".\n";
}


void Particle_Handler::executeParticleBoundaryConditions(){
	for(int i = 0; i < 6; i++){
		// determine whether particles are ghost
        // place ghost particles
        // change the number of particles np_ in each domain
		int inc = boundaries_[i]->computeParticleBCs(&parts_);
		incrementNParticles(inc);
	}
}
#undef _USE_MATH_DEFINES
