#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
//#include "domain.h" 
#include "./IO/IO.hpp"
#include "./grid/grid.hpp"
#include "./particles/particle.hpp"

#if USE_MPI
    #include "mpi.h"  
#else
    #include<time.h>
#endif

int main(int argc, char *argv[]){

    int size,rank;

#if USE_MPI
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    double begin = MPI_Wtime();
    if(rank==0){
        printf("This simulation uses MPI domain decomposition.\n");
    }
#else
    size=1;
    rank=0;
    clock_t begin=clock();
    printf("This simulation is serial.\n");
#endif

    
    // Read command line input //
    if(argc!=2 && rank==0){
      fprintf(stderr,"The correct usage is:\n");
      fprintf(stderr,"  ./EMOOPIC <inputfile>\n");
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
#else
      exit(1);
#endif 
    }

    // Check input file exist 
    if(rank==0){
      FILE *fp = fopen(argv[1],"r");
      if(fp==NULL){
         fprintf(stderr,"Cannot open input file!\n");
#if USE_MPI
         MPI_Abort(MPI_COMM_WORLD,1);
#else
         exit(1);
#endif
      }
      fclose(fp);
    } 


    // Initialization //
    //Domain_t *domains = new Domain(size);
    Grids_list_t *grids; //store Ei,Bi,Rhoi,Ji 
    Particle_Field_list_t *parts_fields; //store x,v,Ex,Bx
    Input_Info_t input_info;
 
    if(rank==0){readinput(argv[1],input_info);}
#if USE_MPI
    Input_Type_t *itype = new_input_type();
    MPI_Datatype infotype; // new type

    MPI_Type_create_struct(itype->count,itype->len,itype->disp,
                           itype->types,&infotype);
    MPI_Type_commit(&infotype);
    MPI_Bcast(&input_info,1,infotype,0,MPI_COMM_WORLD);
#endif

/*
    domain.setup(inputinfo);
    particle.load(inputinfo);//allow restart
    grid.deposeRhoJ(parts);
    grid.poisson(inputinfo); //allow restart
    grid.interpEB(parts);

    // Advance time step //
    t=inputinfo->t0; //initial time
    for(ti=0;ti<nt;ti++){
       particle.dtmin(dt); 
       Pusher.step(part,field,dt);
       particle.pass(domains); //MPI
       grid.deposeRhoJ(parts);
       grid.advanceEB(dt);
       grid.boundary(domains); //MPI
       grid.interpEB(parts);
       // check and write restart files
       if(ti%ntcheck==0){check(t,domains,grids,parts);}
       t+=dt;
     }  

    // output, finalize //
    writeoutput(t,domains,grids,parts); //MPI
*/

#if USE_MPI
    double time = MPI_Wtime()-begin;
    double maxtime;
    MPI_Reduce(&time,&maxtime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
#else
    clock_t end = clock();
    double maxtime = (double)(end-begin)/CLOCKS_PER_SEC;
#endif

    if(rank==0){    
        printf("Program completed successfully!\n");
        printf("Elapsed: %f seconds\n",maxtime);
    }

}
