#include "mpi.h"  
#include "domain.h" 
#include "grid.h"
#include "particle.h"

int main(int argc, char *argv[]){
    MPI_Init(&argc,&argv); //command input: nx,np,nt
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /* Initialization */
    Domain_t *domains = new Domain(size);
    Grid_t *grids = new Grids(nx); //store Ei,Bi,Rhoi,Ji 
    Part_t *parts = new Parts(np); //store x,v,Ex,Bx
    if(rank==0){readinput("input.txt",inputinfo);}
    MPI_Bcast(*inputinfo,1,infotype,0,MPI_COMM_WORLD);

    domain.setup(inputinfo);
    particle.load(inputinfo);//allow restart
    grid.deposeRhoJ(parts);
    grid.poisson(inputinfo); //allow restart
    grid.interpEB(parts);

    /* Advance time step */
    t=inputinfo->t0; //initial time
    for(ti=0;ti<nt;ti++){
       particle.dtmin(dt); 
       particle.push(dt);
       particle.pass(domains); //MPI
       grid.deposeRhoJ(parts);
       grid.advanceEB(dt);
       grid.boundary(domains); //MPI
       grid.interpEB(parts);
       // check and write restart files
       if(ti%ntcheck==0){check(t,domains,grids,parts);}
       t+=dt;
     }  

    /* output, finalize */
    writeoutput(t,domains,grids,parts); //MPI
    MPI_Finalize();
}
