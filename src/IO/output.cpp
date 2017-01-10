#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../particles/particle.hpp"

//#include "../globals.hpp"
#include "output.hpp"

//void writeoutput(Grid *grids, Particle_Handler *parts_fields){
//    //printf("rank %d: writing output files...\n",rank);
//}

//! Check MPI communication by printing data to file
void checkMPI(const char filestem[], double *buffer, int len){
    
    fprintf(stderr,"rank=%d: enter checkMPI\n",rank_MPI);
    assert(buffer!=NULL);

    char fname[50];
    sprintf(fname,"%s%d",filestem,rank_MPI);
    //sprintf(fname,"test%d.txt",rank_MPI);
    fprintf(stderr,"rank=%d: fname=%s\n",rank_MPI,fname);

    FILE *fp;
    fp=fopen(fname,"w+");
    assert(fp!=NULL);
    fprintf(stderr,"rank=%d: file is opened!\n",rank_MPI);
    int i;
    fprintf(stderr,"rank=%d: check looping...\n",rank_MPI);
    for(i=0;i<len;i++){
        //fprintf(stderr,"rank=%d: buffer[%d]=%f\n",rank_MPI,i,buffer[i]);
        fprintf(fp,"%15.8f\n",buffer[i]);
    }
    assert(fclose(fp)==0);
    
    fprintf(stderr,"rank=%d: leaving checkMPI\n",rank_MPI);
}

OutputBoxQuantities::OutputBoxQuantities(Grid* grid, Particle_Handler* handler, Input_Info_t* info) 
	: 	grid_(grid),
		pHandler_(handler)
{
	dStep_		=10;		
	dT_			=-1;
	nextStep_	=0;
	nextT_		=0;
}

void OutputBoxQuantities::output(double t, long i){
	bool needsOutput = false;
	
	if(dT_ > 0 && t >= nextT_){
		needsOutput = true;
		nextT_ += dT_;
	}

	if(dStep_ > 0 && i >= nextStep_){
		needsOutput = true;
		nextStep_ += dStep_;
	}

	if(needsOutput){
		FILE* fp = fopen("history.dat","a");
		double pEn  = 0.0;
		double pM_x = 0.0;			
		double pM_y = 0.0;			
		double pM_z = 0.0;			
		std::vector<Particle> *parts = pHandler_->getParticleVector();
		for(std::vector<Particle>::iterator iter = parts->begin(); iter!=parts->end(); ++iter){
			pEn += 0.5*iter->m *(iter->v[0]*iter->v[0] + iter->v[1]*iter->v[1] + iter->v[2]*iter->v[2]);
			pM_x += iter->m *iter->v[0];
			pM_y += iter->m *iter->v[1];
			pM_z += iter->m *iter->v[2];
		}

		fprintf(fp,"%e %.15e %.15e %.15e %.15e \n",t, pEn,pM_x, pM_y,pM_z);
		fclose(fp);
	}
}

