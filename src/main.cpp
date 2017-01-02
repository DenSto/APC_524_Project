//! Driver program for EMOOPIC
/*! ****************************************************** \n
 * To run the program use the following syntex             \n
 *   Serial version:                                       \n
 *       ./EMOOPIC <inputfile>                             \n
 *   MPI version                                           \n
 *       mpirun -np <nproc> ./EMOOPIC <inputfile>          \n
 *                                                         \n
 * The inputs are specified in <inputfile>                 \n
 **********************************************************/
#define MAIN_CPP
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include <iostream>


#include "./globals.hpp"
#include "./IO/input.hpp"
#include "./IO/output.hpp"
#include "./domain/domain.hpp"
#include "./grid/grid.hpp"
#include "./particles/particle.hpp"
#include "./particles/particle_handler.hpp"
#include "./particles/particle_utils.hpp"
#include "./boundaries/bc_factory.hpp"
#include "./boundaries/boundary_particles.hpp"
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
    size_MPI=size;
    rank_MPI=rank;

    
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

    /* Read and broadcast input file **********************/
    Input *input =  new Input();
    // Master read input file 
    if(rank==0){
      printf("Master reading input file...\n");
      int err = input->readinfo(argv[1]);
      // Check input self-consistency
      err += input->checkinfo();
      if(err!=0) {
        std::cerr << "Input Error. Terminating..." << std::endl;
#if USE_MPI
         MPI_Abort(MPI_COMM_WORLD,1);
#else
         exit(1);
#endif
      }
    }

    Input_Info_t *input_info = input->getinfo();
#if USE_MPI
    // Master broadcast input info
    if(rank==0)printf("Master broadcasting input infomation...\n");
    input->passinfo();
#endif
    debug = input_info->debug; // global debug flag
    int restart = input_info->restart;
    if(debug>1) checkinput(input_info);

    /* Initial setup **************************************/
    if(rank==0)printf("Initial set up...\n");
    // Domain decomposition
    Domain *domain = new Domain(input_info);
    if(debug>1) checkdomain(domain);

    // Initialize particles
    Particle_Handler *part_handler = new Particle_Handler(); 
    part_handler->setPusher(new Boris());

    // Set up particle boundary conditions
    BC_Particle** bc = BC_Factory::getInstance().constructConditions(domain,input_info->parts_bound);
    part_handler->setParticleBoundaries(bc);
    if(debug) fprintf(stderr,"rank=%d:Finish assigning boundary condition\n",rank);

    // Initialize grid
    Grid* grids = new Grid(domain->getnxyz(),domain->getnGhosts(),
               domain->getxyz0(),domain->getLxyz()); //store Ei,Bi,Ji 
    if(debug) fprintf(stderr,"rank=%d: Finish grid constructor\n", rank);

    // Load particles, allow restart
    part_handler->Load(input_info,domain);
    if(debug) fprintf(stderr,"rank=%d: Finish loading particles\n",rank);   

    // Deposite initial charge and current from particles to grid
    //part_handler->depositJ(grids);
    //if(debug) fprintf(stderr,"rank=%d: Finish initial deposition\n",rank);   

    // Solve initial fields from particle or read restart file
    grids->InitializeFields(restart); 
    if(debug) fprintf(stderr,"rank=%d: Finish initializing fields\n",rank);   

    // Interpolate fields from grid to particle
    // Prepare initial push of particles
    part_handler->InterpolateEB(grids);
    if(debug) fprintf(stderr,"rank=%d: Finish initializing interpolation\n",rank);   

    // prepare ghost cells: either MPI neighbors or physical boundary 
    int xgsize = grids->getGhostVecSize();
    int ygsize = 1; //dummy
    int zgsize = 1; //dummy
    domain->mallocGhosts(xgsize,ygsize,zgsize);
    if(debug) fprintf(stderr,"rank=%d: Finish allocating ghosts\n",rank);   

    // prepare time step
    int nt = input_info->nt; //number of steps to run
    time_phys = input_info->t0; //initial time
    double dt = 1/domain->getmindx(); //c=1, resolve EM wave
    if(debug) fprintf(stderr,"rank=%d: Finish preparing time step\n",rank);   

    /* Advance time step **********************************/
    if(rank==0)printf("Advancing time steps...\n");
    for(int ti=0;ti<nt;ti++){
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Time Loop\n",rank,ti);   
       // check and write restart files
//       if(ti%ntcheck==0){check(t,domains,grids,parts);}

       // push particles
       part_handler->Push(dt);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Push\n",rank,ti);   

       // Pass particle through MPI boundary, or physical boundary conditions
//       part_handler->executeParticleBoundaryConditions();
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass parts\n",rank,ti);   

       // deposite charge and current on grid
       part_handler->depositRhoJ(grids,false);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish deposition\n",rank,ti);   

       // evolve E, B fields
       grids->evolveFields(dt);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish evolve\n",rank,ti);   

       // pass field boundaries 
       domain->PassFields(grids,input_info,-1);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass fields\n",rank,ti);   

       // Interpolate fields from grid to particle
       part_handler->InterpolateEB(grids);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish interpolate\n",rank,ti);   

       // remove any particles left in the ghost cells
       part_handler->clearGhosts();

       time_phys += dt;
     }  


    /* output, finalize ***********************************/
    if(rank==0)printf("Writing output files...\n");
    writeoutput(grids,part_handler); //MPI
    if(debug) fprintf(stderr,"rank=%d: Finish writeoutput\n",rank);   

    // free memory
    domain->freeGhosts(); // Ghost for field MPI
    delete domain; 
    delete [] bc; // particle boundary condition
    delete part_handler;
    delete grids;
    delete input;
    if(debug) fprintf(stderr,"rank=%d: Finish free\n",rank);   

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
