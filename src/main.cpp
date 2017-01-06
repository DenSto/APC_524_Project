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
#include "./domain/domain.hpp"
#include "./grid/grid.hpp"
#include "./poisson/poisson.hpp"

#include "./IO/input.hpp"
#include "./IO/output.hpp"

#include "./particles/particle.hpp"
#include "./particles/particle_handler.hpp"
#include "./particles/particle_utils.hpp"

#include "./boundaries/particle_bc_factory.hpp"
#include "./boundaries/particles_boundary.hpp"
#include "./boundaries/field_bc_factory.hpp"
#include "./boundaries/fields_boundary.hpp"

#include "./pusher/pusher.hpp"
#include "./pusher/boris.hpp"
#include "./pusher/relativisticBoris.hpp"

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
      // Check input self-consistency and load physical units to input
      err += input->checkinfo(); // should not be commented out
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
    int restart = input_info->restart;
    int relativity = input_info->relativity;
    debug = input_info->debug; // global debug flag
    if(debug>1) checkinput(input_info);

    /* Initial setup **************************************/
    if(rank==0)printf("Initial set up...\n");
    // Domain decomposition
    Domain *domain = new Domain(input_info);
    if(debug>1) checkdomain(domain);

    // Initialize particles and pusher
    Particle_Handler *part_handler = new Particle_Handler(); 
    if(relativity==0){
       if(rank==0)printf("    Use non-relativistic pusher\n");
       part_handler->setPusher(new Boris());
    }else{
       if(rank==0)printf("    Use relativistic pusher\n");
       part_handler->setPusher(new Relativistic_Boris());
    } 

    // Set up particle boundary conditions
    BC_Particle** bc = Part_BC_Factory::getInstance().constructConditions(domain,input_info->parts_bound);
    part_handler->setParticleBoundaries(bc);
    if(debug) fprintf(stderr,"rank=%d:Finish assigning particle boundary condition\n",rank);

    // Initialize grid
    Grid *grids;
    if(restart>0){//no need to solve Poisson's equation
        if(rank==0)printf("    Grid initialing...\n");
        grids = new Grid(domain->getnxyz(),domain->getnGhosts(),
               domain->getxyz0(),domain->getLxyz()); //store Ei,Bi,Ji 
    }else{//need to solve Poisson's equation
        if(rank==0)printf("    Grid initialing: will solve Poisson's equations...\n");
        grids = new Poisson_Solver(domain,input_info);
    }
    if(debug) fprintf(stderr,"rank=%d: Finish grid constructor\n", rank);

    // Set up field boundary conditions
    Field_BC_Factory::getInstance().Construct(domain,grids,input_info);
    if(debug) fprintf(stderr,"rank=%d:Finish assigning field boundary condition\n",rank);

    // Load particles, allow restart
    if(rank==0)printf("    Loading particles...\n");
    part_handler->Load(input_info,domain);
    if(debug) fprintf(stderr,"rank=%d: Finish loading particles\n",rank);   

    // if initial run, Deposite charge and current from particles to grid
    if(restart==0){
        if(rank==0)printf("    Depositing rho and J for Poisson solver...\n");
        part_handler->depositRhoJ(grids,true,domain,input_info);
        if(debug) fprintf(stderr,"rank=%d: Finish initial deposition\n",rank);   
    }

    // Solve initial fields from particle or read restart file
    if(rank==0)printf("    Initializing fields...\n");
    grids->InitializeFields(); 
    if(debug) fprintf(stderr,"rank=%d: Finish initializing fields\n",rank);   

    // Interpolate fields from grid to particle
    // Prepare initial push of particles
    part_handler->InterpolateEB(grids);
    if(debug) fprintf(stderr,"rank=%d: Finish initializing interpolation\n",rank);   

    // prepare ghost cells: either MPI neighbors or physical boundary 
    int xgsize = grids->getGhostVecSize(-1);
    int ygsize = 1; //dummy
    int zgsize = 1; //dummy
    domain->mallocGhosts(xgsize,ygsize,zgsize);
    if(debug) fprintf(stderr,"rank=%d: Finish allocating ghosts\n",rank);   

    // prepare time step
    int nt = input_info->nt; //number of steps to run
    time_phys = input_info->t0; //initial time
    dt_phys = 1/domain->getmindx(); //c=1, resolve EM wave
    if(debug) fprintf(stderr,"rank=%d: Finish preparing time step\n",rank);   

    // write initial restart files
    // write initial diagnostic files
    part_handler->outputParticles(0,input_info);

    /* Advance time step **********************************/
    if(rank==0)printf("Advancing time steps...\n");
    for(int ti=0;ti<nt;ti++){
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Time Loop\n",rank,ti);   
       // push particles
       part_handler->Push(dt_phys);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Push\n",rank,ti);   

       // Pass particle through MPI boundary, or physical boundary conditions
//       part_handler->executeParticleBoundaryConditions();
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass parts\n",rank,ti);   

       // deposite charge and current on grid
//       part_handler->depositRhoJ(grids,false,domain,input_info);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish deposition\n",rank,ti);   

       // evolve E, B fields
       grids->evolveFields(dt_phys);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish evolve\n",rank,ti);   

       // pass field boundaries 
       grids->executeBC();
       //domain->PassFields(grids,input_info,-1);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass fields\n",rank,ti);   

       // Interpolate fields from grid to particle
       part_handler->InterpolateEB(grids);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish interpolate\n",rank,ti);   

       // remove any particles left in the ghost cells
       part_handler->clearGhosts();
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish clearGhosts\n",rank,ti);   

       time_phys += dt_phys;

       // check and write restart files
//       if(ti%ntcheck==0){check(t,domains,grids,parts);}
       // output diagnostic files
       part_handler->outputParticles(ti+1,input_info);
     }  

    //Output particle velocities
    //part_handler->outputParticleVel();

    /* output, finalize ***********************************/
    if(rank==0)printf("Writing output files...\n");
    writeoutput(grids,part_handler); //MPI
    if(debug) fprintf(stderr,"rank=%d: Finish writeoutput\n",rank);   

    // free memory
    domain->freeGhosts(); // Ghost for field MPI
    delete domain; 
    delete [] bc; // particle boundary condition
    delete part_handler;
    grids->freeBoundaries(); // field boundary conditions
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
