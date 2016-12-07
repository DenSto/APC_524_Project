//! Driver program for EMOOPIC
/*! ******************************************************
 * To run the program use the following syntex
 *   Serial version:
 *       ./EMOOPIC <inputfile>
 *   MPI version
 *       mpirun -np <nproc> ./EMOOPIC <inputfile>
 *
 * The inputs are specified in <inputfile>
 **********************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

#include "./IO/IO.hpp"
#include "./domain/domain.hpp"
#include "./grid/grid.hpp"
#include "./particles/particle.hpp"
#include "./particles/particle_list.hpp"
#include "./particles/particle_utils.hpp"
#include "./pusher/pusher.hpp"
#include "./pusher/boris.hpp"

#if USE_MPI
    #include "mpi.h"  
#else
    #include<time.h>
#endif

int main(int argc, char *argv[]){

    int size,rank=0;

    /* Initialize *****************************************/
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

    
    /* Read and check command line input ******************/
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


    /* Read and broadcast input file **********************/
    Input_Info_t input_info;
    if(rank==0){readinput(argv[1],&input_info);}
#if USE_MPI
    Input_Type itype;
    MPI_Datatype infotype; // new type

    MPI_Type_create_struct(itype.getcount(),itype.getlens(),itype.getdisps(),
                           itype.gettypes(),&infotype);
    MPI_Type_commit(&infotype);
    MPI_Bcast(&input_info,1,infotype,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    checkinput(rank,&input_info);

    /* Initial setup **************************************/
    // Domain decomposition
    Domain *domain = new Domain(size,rank,&input_info);
    checkdomain(rank,domain);
    //domain.setup(inputinfo);

    // Initialize particles
    Particle_Field_List *parts_fields = new Particle_Field_List(input_info.np); 
    parts_fields->setPusher(new Boris());

    // Initialize grid
    Grid* grids = new Grid(domain->getnxyz(),domain->getnGhosts(),
               domain->getxyz0(),domain->getLxyz()); //store Ei,Bi,Ji 

    // Load particles, allow restart
    parts_fields->Load();

    // Deposite initial charge and current from particles to grid
    parts_fields->depositCharge(grids);
    parts_fields->depositCurrent(grids);

    // Solve initial fields from particle or read restart file
    grids->InitializeFields(&input_info); 

    // Prepare to push particles
    parts_fields->InterpolateEB(grids);

    /* Advance time step **********************************/
    double t = input_info.t0; //initial time
/*    for(ti=0;ti<nt;ti++){
       particle.dtmin(dt); 
       Pusher.step(part,field,dt);
       particle.pass(domains); //MPI
       grid.depositRhoJ(parts_fields);
       grid.advanceEB(dt);
       grid.boundary(domains); //MPI
*/
       parts_fields->InterpolateEB(grids);
/*
       // check and write restart files
       if(ti%ntcheck==0){check(t,domains,grids,parts);}
       t+=dt;
     }  

    // output, finalize //
    writeoutput(t,domains,grids,parts); //MPI
*/
    delete domain;
    delete parts_fields;
    delete grids;

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
#if USE_MPI
	MPI_Finalize();
#endif

}
