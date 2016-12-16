#include <stdio.h>
#include <assert.h>
#include "../globals.hpp"
#include "output.hpp"

void writeoutput(double t, int rank, Grid *grids, Particle_Handler *parts_fields){
    //printf("rank %d: writing output files...\n",rank);
}

//! Check MPI communication by printing data to file
void checkMPI(const char filestem[], double *buffer, int len){
    
    fprintf(stderr,"rank=%d: enter checkMPI\n",rank_MPI);
    assert(buffer!=NULL);

    char fname[50];
    strcpy(fname,filestem);
    char append[20];
    sprintf(append,"%d",rank_MPI);
    fprintf(stderr,"rank=%d: fname=%s,append=%s\n",rank_MPI,fname,append);
    strcat(fname,append);
    fprintf(stderr,"rank=%d: fname=%s\n",rank_MPI,fname);

    FILE *fp;
    fp=fopen(fname,"w+");
    assert(fp!=NULL);
    int i;
    fprintf(stderr,"rank=%d: check looping...\n",rank_MPI);
    for(i=0;i<len;i++){
        //fprintf(stderr,"rank=%d: buffer[%d]=%f\n",rank_MPI,i,buffer[i]);
        fprintf(fp,"%15.8f\n",buffer[i]);
    }
    fclose(fp);
    
    fprintf(stderr,"rank=%d: leaving checkMPI\n",rank_MPI);
}

