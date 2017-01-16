#include "../globals.hpp"
#include "particle_handler.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"
#include "../grid/grid.hpp"
#include "../utils/RNG.hpp"
#include <string.h>
#include <iostream>
#include <limits>
#define _USE_MATH_DEFINES
#if USE_MPI
#include "mpi.h"
#endif

Particle_Handler::Particle_Handler(){
	np_=0;
	dT_=-1;
	dstep_=-1;
	nextT_=0;
	nextStep_=0;
}

Particle_Handler::~Particle_Handler(){
}

//! Load and initialize the particle handler. Should be called at the beginning of the run.
void Particle_Handler::Load(Input_Info_t *input_info, Domain* domain){

    int restart = input_info->restart;
    double* L = domain->getLxyz();
    double* x0 = domain->getxyz0();

    int nspec      = input_info->nspecies;
    int npart      = input_info->nparticles_domain;
    double *mass   = input_info->mass_ratio;
    double *charge = input_info->charge_ratio;
    double *dens   = input_info->dens_frac; 
    int *test   = input_info->isTestParticle; 

    if(restart==0){// initial run
		if(rank_MPI==0)printf("    Loading random particles from distribution...\n");
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
			assert(ispec<nspec);
			p.q = charge[ispec]; // super particle charge
			p.m = mass[ispec];   // super particle mass
			p.isTestParticle = test[ispec];
			p.type = ispec;

            if(debug>3)fprintf(stderr,"charge=%f\n",p.q);

			vth=UNIT_VTH*sqrt(input_info->temp[ispec]/p.m);

			p.x[0]=rng->getUniform()*L[0]+x0[0];
			p.x[1]=rng->getUniform()*L[1]+x0[1];
			p.x[2]=rng->getUniform()*L[2]+x0[2];

			p.v[0]=rng->getGaussian(0.0,vth);
			p.v[1]=rng->getGaussian(0.0,vth);
			p.v[2]=rng->getGaussian(0.0,vth);

			p.gamma = sqrt(1 - p.v[0]*p.v[0] - p.v[1]*p.v[1] - p.v[2]*p.v[2]);

			p.my_id=ip;
			p.initRank=rank_MPI;

			parts_.push_back(p);
			np_++;
    	} 
	}
    else{//read restart file
	if(rank_MPI==0)printf("    Loading particles from file...\n");
        //dummy code inserted by YShi for testing
        //insert a single particle at the center of the cell
	Particle p = new_particle();
	p.q = 1.0;
	p.m = 1.0;
	p.x[0]=L[0]/2+x0[0];
	p.x[1]=L[1]/2+x0[1];
	p.x[2]=L[2]/2+x0[2];
	p.v[0]=0.01;
	p.v[1]=0.001;
	p.v[2]=0.001;
	parts_.push_back(p);
	np_++;
    }

    //if boundaries are ALL periodic, and initialization is Poisson, then charge and velocity must be zero (within each species)
    bool allBoundariesPeriodic = true;
    for (int i=0; i<6; i++) {
      if ( strcmp(input_info->fields_bound[i],"periodic") != 0 ) {
	allBoundariesPeriodic = false;
	break;
      }
    }
    if (allBoundariesPeriodic && strcmp(input_info->fields_init,"poisson") == 0 ) {
      //zero them out.
      printf("============================================================================================\n");
      printf("NOTICE: Because the user requested Poisson initialization with periodic boundary conditions,\n");
      printf("        the average particle velocity will be initialized to zero (by species).\n");
      printf("============================================================================================\n");
      zeroAvgChargeAndVelocity_();
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
  Interpolator *interpolator = new Interpolator();

  long iCell = -1; //cell # tracker
  long pCell = 0; //particle cell #
  double cellvars[21];//Vector describing position of and all field elements of a cell
                      //["least" corner vertex, E-field on edges, B-field on surfaces]
                      //21 elements ordered as: [x[0],x[1],x[2],E1,E2,E3,B1,B2,B3]
                      //              of sizes:   1, 1, 1, 4, 4, 4, 2, 2, 2
  double pos[3]; //Vector of position of particle.
  double lcell[3]; //Vector of lengths of unit cell.

  //Get lengths of grid cells.
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  for (long i=0; i<np_; i++) {
	assert(!parts_[i].isGhost);
	
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
	if(debug>3)  fprintf(stderr,"rank=%d,call Interpolator\n",rank_MPI);
    interpolator->interpolate_fields(pos, lcell, cellvars, &(parts_[i].field));
  }
  if(debug>3) fprintf(stderr,"rank=%d,Finish InterpolateEB\n",rank_MPI);
}

//! Sort particles based on grid location. 
/*
 * Sorts particles based on grid location using std::sort which is O(n log n).
 * This should be called often to ensure cache hits.  Should be emperically determined.
 */
void Particle_Handler::SortParticles(Particle_Compare comp){
	std::sort(parts_.begin(),parts_.end(),comp);
}

void Particle_Handler::depositRhoJ(Grid *grid, bool depositRho, Domain* domain, Input_Info_t* input_info){
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
  grid->constJ(0,0,0);
  if (depositRho) grid->constRho(0);

  //Cycle through particles, depositing RhoJ for each one.
  for (long i=0; i<np_; i++) {
    //Only if particle ENDS in a 'real' (non-ghost) cell and is not a test particle, do we 'deposit' it.
    if (!parts_[i].isGhost && !parts_[i].isTestParticle){
    	//Get particle's 'entered' (post-push) cell id
		tempCellID = grid->getCellID(parts_[i].x[0],parts_[i].x[1],parts_[i].x[2]);

	    //If list of particles continues to be in same cell, add current to existing JObj object.
    	//Otherwise, deposit (add) RhoJ to the grid, and re-point (and zero out) the existing JObj object.
		if (tempCellID != cellID) {
			//If cellID has already been assigned...
			if (cellID != -1) {
				//Deposit to grid
				grid->addJ(cellID,JObj);
				if (depositRho) grid->addRho(cellID,RhoObj);
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
  grid->addJ(cellID,JObj);
  if (depositRho) grid->addRho(cellID,RhoObj);

  //main will pass rho and J fields between domains.
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

//! Clear all ghost particles. Uses a swap-to-back and pop-last-element for speed.
void Particle_Handler::clearGhosts(){
        int nghost = 0;
	for(std::vector<Particle>::iterator iter = parts_.begin(); iter != parts_.end();){
		if(iter->isGhost){
			std::swap(*iter, parts_.back());
			parts_.pop_back();
                        nghost +=1;
		} else {
			iter++;
		}
	}
	assert((long)parts_.size() == np_);
        if(debug>2)fprintf(stderr,"rank=%d: %d ghosts are cleared.\n",rank_MPI,nghost);
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

//! Output particles
/*!
 *	Output particles. Currently outputs time, position and velocity.
 *
 *	Can work with either cadencing on time (output every dT) or 
 *	cadencing on steps (output ever dsteps), or both. Either are 
 *	optional parameters in the input file and will default to -1.
 *
 *  Particles are written to the same directory as the program input file and 
 *  named with initial rank and id.
 */
void Particle_Handler::outputParticles(const char* basename,long step, Input_Info_t *input_info){

	if(debug>2 && rank_MPI==0) printf("    ti=%ld: writing particle tracks...\n",step);

    double t = time_phys;
    dstep_ = input_info->nstep_parts;
    outputCount_ = input_info->output_pCount;

	static bool init = true;
	bool needsOutput=false;
	if(init){
		init=false;
		//mkdir("./tracks", 0775); //create the particle directory
		for(std::vector<Particle>::iterator iter = parts_.begin();iter!=parts_.end();++iter){
			if(iter->my_id < outputCount_){
				char fname[100];
				sprintf(fname,"%strack_%d_%ld.dat",basename,iter->initRank,iter->my_id);	
				FILE *pout=fopen(fname,"w");
				fprintf(pout,"[1] time [2] x [3] y [4] z  [5] vx   [6] vy   [7] vz\n");
				fclose(pout);
			}
		}
	}

	// cadence on i or t
	if(dT_ > 0 && t  >= nextT_){
		needsOutput=true;
		nextT_ += dT_;
	}
	
	if(dstep_ > 0 && step >= nextStep_){
		needsOutput=true;
		nextStep_+=dstep_;
	}

	if(!needsOutput)
		return;

        if(debug){
            fprintf(stderr,"rank=%d: writing tracks for %ld particles...\n",
                            rank_MPI,outputCount_);
        }
	char fname[100];
	FILE *pout;
	for(std::vector<Particle>::iterator iter = parts_.begin();iter!=parts_.end();++iter){
		if(iter->my_id < outputCount_){
			sprintf(fname,"%strack_%d_%ld.dat",basename,iter->initRank,iter->my_id);	
                        if(debug>1)fprintf(stderr,"    track file name %s\n",fname);
			pout=fopen(fname,"a");
			assert(pout != NULL);
			fprintf(pout,"%e %.15e %.15e %.15e %.15e %.15e %.15e\n",t,
					iter->x[0],iter->x[1],iter->x[2],
					iter->v[0],iter->v[1],iter->v[2]);
			fclose(pout);
		}
	}
        if(debug)fprintf(stderr,"rank=%d: finish writing particle tracks!\n",rank_MPI);
}

void Particle_Handler::outputParticleVel(){
  for (long i=0; i<np_; i++) {
    //Get velocity of particle and species
    char fname[100];
    FILE *fp;
    sprintf(fname,"./velocity_%d.dat",rank_MPI);
    fp = fopen(fname,"w+");
    fprintf(fp,"Species = %f   Vel = %f\n",parts_[i].q,parts_[i].v[0]);
    fclose(fp);
  }
}

void Particle_Handler::zeroAvgChargeAndVelocity_(){

  double qTot = 0.0;   //total charge of all particles
  long qRefIdx = 0;    //reference index
  long qCount = 0;     //count of the # of charges
  double vx = 0.0;     //x velocity
  double vy = 0.0;     //y velocity
  double vz = 0.0;     //z velocity

  double min_q = std::numeric_limits<float>::max();  //minimum absolute charge on any particle
  std::vector<int> countedParticle(np_,0);           //vector of integers, tracking the handling of particles

  do {
    qRefIdx = 0;
    qCount = 0;
    vx = 0.0;
    vy = 0.0;
    vz = 0.0;

    //Sum up velocities of particles with charge matching qRef.
    for (long i=0; i<np_; i++) {
      if (countedParticle[i]==0) {
	if (qRefIdx == 0) qRefIdx = i;

	if (parts_[i].q == parts_[qRefIdx].q) {
	  countedParticle[i] = 1; //move to 'state 1: counted velocity toward an average'
	  qCount += 1;
	  vx += parts_[i].v[0];
	  vy += parts_[i].v[1];
	  vz += parts_[i].v[2];
	} 
      }
    }

    if (qCount > 0) {
      //Take the average
      vx /= qCount;
      vy /= qCount;
      vz /= qCount;

      //Subtract out the average
      for (long i=0; i<np_; i++) {
	if (countedParticle[i]==1) {
	  countedParticle[i] = 2; //move to 'state 2: average deducted from velocity, and added charge to qTot'
	  parts_[i].v[0] -= vx;
	  parts_[i].v[1] -= vy;
	  parts_[i].v[2] -= vz;
	  qTot += parts_[i].q;
	  min_q = std::min(fabs(parts_[i].q), min_q);
	}
      }
    }

  } while(qCount > 0);

  //Test for errors
  if(fabs(qTot) >= min_q / np_) {
    printf("ERROR: Poisson initialization with periodic boundaries cannot be run with a net charge.  Please adjust inputs and rerun.\n");
  }
  assert(fabs(qTot) < min_q / np_);
  for (long i=0; i<np_; i++) {
    assert(countedParticle[i] == 2);
  }

  countedParticle.clear();
}

#undef _USE_MATH_DEFINES
