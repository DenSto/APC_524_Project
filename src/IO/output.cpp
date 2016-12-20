#include <stdio.h>
#include <assert.h>

#include "../globals.hpp"
#include "output.hpp"

void writeoutput(Grid *grids, Particle_Handler *parts_fields){
    //printf("rank %d: writing output files...\n",rank);
}

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

