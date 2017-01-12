#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../particles/particle.hpp"

#include "../globals.hpp"
#include "output.hpp"

#define SQR(x) (x*x)

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

OutputBoxQuantities::OutputBoxQuantities(const char* fname,Grid* grid, Particle_Handler* handler, Input_Info_t* info) 
	: 	grid_(grid),
		pHandler_(handler)
{
	sprintf(filename_,"%s",fname);
	dStep_		=10;		
	dT_			=-1;
	nextStep_	=0;
	nextT_		=0;
}

void OutputBoxQuantities::output(double t, long i){
	bool needsOutput = false;
	static bool init = true;

	if(init){
		init = false;
		FILE* fp = fopen(filename_,"w");
		fprintf(fp,"[1] time [2] part mom x [3] part mom y [4] part mom z [5] part Energy [6] B Energy [7] E Energy [8] Total E\n");
		fclose(fp);
	}
	
	if(dT_ > 0 && t >= nextT_){
		needsOutput = true;
		nextT_ += dT_;
	}

	if(dStep_ > 0 && i >= nextStep_){
		needsOutput = true;
		nextStep_ += dStep_;
	}

	if(needsOutput){
		FILE* fp = fopen(filename_,"a");

		// Particle energy and momentum
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

		// field energy
		int beg[6];
		double EEx = 0;
		double EEy = 0;
		double EEz = 0;
		double EBx = 0;
		double EBy = 0;
		double EBz = 0;

		double vol = grid_->getCellVolume();
		double**** fields = grid_->getFieldPtr();

		grid_->getRealIndices(0,beg);
		for(int i = beg[0]; i < beg[3]; i++){
			for(int j = beg[1]; j < beg[4]; j++){
				for(int k = beg[2]; k < beg[5]; k++){
					EEx += 0.5*vol*SQR(0.5*(fields[E_X][i][j][k] + fields[E_X][i+1][j][k]));
					EEy += 0.5*vol*SQR(0.5*(fields[E_Y][i][j][k] + fields[E_Y][i][j+1][k]));
					EEz += 0.5*vol*SQR(0.5*(fields[E_Z][i][j][k] + fields[E_Z][i][j][k+1]));
					EBx += 0.5*vol*SQR(0.5*(fields[B_X][i][j][k] + fields[B_X][i+1][j][k]));
					EBy += 0.5*vol*SQR(0.5*(fields[B_Y][i][j][k] + fields[B_Y][i][j+1][k]));
					EBz += 0.5*vol*SQR(0.5*(fields[B_Z][i][j][k] + fields[B_Z][i][j][k+1]));
				}
			}
		}
		
		double B_en = EBx + EBy + EBz;
		double E_en = EEx + EEy + EEz;
		double E_tot = B_en + E_en + pEn;

		fprintf(fp,"%e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",t, pM_x, pM_y,pM_z,pEn, B_en, E_en, E_tot);
		fclose(fp);
	}
}

#undef SQR
