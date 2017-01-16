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
#include "./domain/resolve.hpp"
#include "./IO/input.hpp"
#include "./IO/output.hpp"
#include "./IO/hdf5io.hpp"
#include "./grid/grid.hpp"
#include "./poisson/poisson.hpp"


#include "./particles/particle.hpp"
#include "./particles/particle_handler.hpp"
#include "./particles/particle_utils.hpp"

#include "./boundaries/particle_bc_factory.hpp"
#include "./boundaries/field_bc_factory.hpp"

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
      err += input->ProcessInfo(); // should not be commented out
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
    int electrostatic = input_info->electrostatic;
    debug = input_info->debug; // global debug flag
    if(debug>1) checkinput(input_info);

    /***************************************************************************/
    /* Initial setup                                                           */
    /***************************************************************************/
    if(rank==0)printf("Initial set up...\n");
    // Domain decomposition
    Domain *domain = new Domain(input_info);
    if(debug>1) checkdomain(domain);

    // determine temporal and spatial resolution
    Resolution *resolution = new Resolution(input_info);

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
    BC_Particle** bc = Part_BC_Factory::getInstance().constructConditions(domain,input_info);
    part_handler->setParticleBoundaries(bc);
    if(debug) fprintf(stderr,"rank=%d:Finish assigning particle boundary condition\n",rank);

    // Initialize grid
    Grid *grids;
    if(restart==0 && strcmp(input_info->fields_init,"poisson")==0){
        //need to solve Poisson's equation
        if(rank==0)printf("    Grid initialing: will solve Poisson's equations...\n");
        grids = new Poisson_Solver(domain,input_info);
    }else{
        //no need to solve Poisson's equation
        if(rank==0)printf("    Grid initialing: won't solve Poisson's equations...\n");
        grids = new Grid(domain->getnxyz(),1,domain->getxyz0(),domain->getLxyz()); 
    }
    if(debug) fprintf(stderr,"rank=%d: Finish grid constructor\n", rank);

    // Set up field boundary conditions
    Field_BC_Factory::getInstance().Construct(domain,grids,input_info);
    if(debug) fprintf(stderr,"rank=%d:Finish assigning field boundary condition\n",rank);

    // Load particles, allow restart
    if(rank==0)printf("    Loading particles...\n");
    part_handler->Load(input_info,domain);
    if(debug) fprintf(stderr,"rank=%d: Finish loading particles\n",rank);   

    // Deposit charge and current from particles to grid
    if(restart==0 && strcmp(input_info->fields_init,"poisson")==0){
        if(rank==0)printf("    Depositing rho and J for Poisson solver...\n");
        part_handler->depositRhoJ(grids,true,domain,input_info);
        // sum charge and current on MPI boundaries 
        grids->executeBC(-2); // boundary for Rho,J
#if USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if(debug) fprintf(stderr,"rank=%d: Finish Pass initial R,J\n",rank);
        if(debug) fprintf(stderr,"rank=%d: Finish initial deposition\n",rank);   
    }

    // Initialize fields by solving poisson or read files
    if(rank==0)printf("    Initializing fields...\n");
    grids->InitializeFields(input_info);
    if(debug) fprintf(stderr,"rank=%d: Finish initializing fields\n",rank);

    // Interpolate fields from grid to particle
    // Prepare initial push of particles
    part_handler->InterpolateEB(grids);
    if(debug) fprintf(stderr,"rank=%d: Finish initializing interpolation\n",rank);

    // initialize outputs
    // format file names
    std::string inputname = argv[1];
    std::string dir = inputname.substr(0,inputname.find_last_of('/')+1);
    //if(dir == "") dir = ".";
    std::string outputname = dir + "output.h5";
    std::string boxName = dir + "history.dat";
    int nstep_fields  = input_info->nstep_fields;
//    int nstep_restart = input_info->nstep_restart;
    int output_fields = input_info->which_fields;
    if(rank==0)printf("    Writing initial field diagnostic files...\n");

    // Box average quantities
    OutputBoxQuantities* boxOutput = new OutputBoxQuantities(boxName.c_str(),grids,part_handler,input_info);
    if(debug) fprintf(stderr,"rank=%d: Finish initialize particle output\n",rank);

    // particle output
    part_handler->outputParticles(dir.c_str(),0,input_info); 
    if(debug) fprintf(stderr,"rank=%d: Finish writing initial particle output\n",rank);

    // fields output
    Hdf5IO* hdf5io = new Hdf5IO(outputname.c_str(),grids,domain,output_fields);
    if(debug) fprintf(stderr,"rank=%d: Finish initialize field output\n",rank);
    if(output_fields>=0){
        hdf5io->writeFields(grids, time_phys);
        if(debug) fprintf(stderr,"rank=%d: Finish writing initial fields output\n",rank);
    }

    // prepare time step
    int nt = input_info->nt; //number of steps to run
    time_phys = input_info->t0; //initial time
    dt_phys = domain->getmindx()/1; //c=1, resolve EM wave
    dt_phys /= 100.0;
    if(debug) fprintf(stderr,"rank=%d: Finish preparing time step\n",rank);

    /***************************************************************************/
    /* Advance time step                                                       */
    /***************************************************************************/
    if(rank==0)printf("Advancing time steps with dt = %f ps\n",dt_phys*UNIT_TIME);
    for(int ti=0;ti<nt;ti++){

       if(rank==0 && ti%100==0)fprintf(stderr,"ti=%d\t\t t=%f ps\n",ti,time_phys*UNIT_TIME);   
		
       boxOutput->output(time_phys,ti);

       /* push particles ***********************/
       part_handler->Push(dt_phys);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Push\n",rank,ti);

       // Pass particle across MPI boundaries, or implement physical boundary conditions
       // All particles are in physical cells, no particle lives in ghost cell
       part_handler->executeParticleBoundaryConditions();
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass parts\n",rank,ti);

       // remove any particles left in the ghost cells
       part_handler->clearGhosts();
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish clearGhosts\n",rank,ti);

       /* deposit charge and current ***********/
       // only deposit particles in physical cells
       part_handler->depositRhoJ(grids,true,domain,input_info);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish deposition\n",rank,ti);   

       // sum charge and current on MPI boundaries 
       grids->executeBC(-2); // boundary for Rho,J
#if USE_MPI
       MPI_Barrier(MPI_COMM_WORLD);
#endif
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass R,J\n",rank,ti);

       /* evolve E, B fields *******************/
       if(electrostatic==1){ // electrostatic field evolve
           grids->evolveFieldsES(dt_phys);
       }else { // electromagnetic field evolve
           grids->evolveFields(dt_phys);
       }
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish evolve\n",rank,ti);

       // pass E,B field across MPI boundaries, or implement physical boundary conditions
       grids->executeBC(-1); //E,B load physical to replace ghost
#if USE_MPI
       MPI_Barrier(MPI_COMM_WORLD);
#endif
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish Pass E,B\n",rank,ti);

       /* interpolation ************************/
       // Interpolate fields from grid to particle
       part_handler->InterpolateEB(grids);
       if(debug>1) fprintf(stderr,"rank=%d,ti=%d: Finish interpolate\n",rank,ti);

       time_phys += dt_phys;

       /* output *******************************/
       // writing files
       if((ti+1)%nstep_fields==0) {
         if(debug && rank==0)printf("    ti=%d: Writing field diagnostic files...\n",ti+1);
         // fields output
         if(output_fields>=0) hdf5io->writeFields(grids, time_phys);
       }
       // particle output
       part_handler->outputParticles(dir.c_str(),ti+1,input_info); 

     } // timestep loop

     if(rank==0) fprintf(stderr,"ti=%d\t\t t=%f ps\n",nt,time_phys*UNIT_TIME);   
     if(rank==0) printf("***Timestep loop complete***\n");

    /***************************************************************************/
    /* output, finalize                                                        */
    /***************************************************************************/
    if(rank==0)printf("Writing final output files...\n");
    if((nt)%nstep_fields!=0) { // if we didn't write on final timestep, write
      if(rank==0)printf("    Writing final field diagnostic files...\n");
      // fields output
      if(output_fields>=0) hdf5io->writeFields(grids, time_phys);
    }
    //writeoutput(grids,part_handler); //MPI
    if(debug) fprintf(stderr,"rank=%d: Finish final output\n",rank);

    //Output particle velocities
    //part_handler->outputParticleVel();

    // free memory
    delete hdf5io;
    delete domain;
    delete [] bc; // particle boundary condition
    delete part_handler;
    grids->freeBoundaries(); // field boundary conditions
    delete grids;
    delete resolution;
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
