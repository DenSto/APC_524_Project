#include <stdio.h>
#include <assert.h>
#include "domain.hpp"
#include "../IO/input.hpp"
#include "../IO/output.hpp"
#include "../globals.hpp"

#if USE_MPI
#include "mpi.h"
#endif

//! Pass fields across MPI boundaries, or execute physical boundary conditions
void Domain::PassFields(Grid *grids, Input_Info_t *input_info, int sendID){
#if USE_MPI
    int tagl=1; // MPI tag: to left, or from right
    int tagr=2; // MPI tag: to right, or from left
#endif
    char (*bound)[NCHAR] = input_info->fields_bound;
    
    /* x field boundaries *********************************/
    int xgsize = grids->getGhostVecSize(); // ghost size in x direction
    assert(xgsize>0);   
    if(debug>1) fprintf(stderr,"rank=%d: xgsize=%d\n",rank_,xgsize);
    int offl=0;
    int offr=xgsize;
    
    // load buffer xghost_send_
    grids->getGhostVec(-1, &(xghost_send_[offl]),sendID);// left
    grids->getGhostVec(1, &(xghost_send_[offr]),sendID);// right
    if(debug>2){
       fprintf(stderr,"rank=%d:checking send\n",rank_);
       checkMPI("xsend.dat",xghost_send_,2*xgsize);
       fprintf(stderr,"rank=%d:finished checking send\n",rank_);
    }

    // left boundary
    if(rank_>rank_xl_ || strcmp(bound[0],"periodic")==0){// left boundary is MPI
       if(debug) fprintf(stderr,"rank=%d: xl is MPI\n",rank_);
#if USE_MPI
       // send to left neighbor
       MPI_Isend(&(xghost_send_[offl]),xgsize,MPI_DOUBLE,rank_xl_,
                 tagl,MPI_COMM_WORLD,&(xreqs_[0]));
       // recv from left neighbor
       MPI_Irecv(&(xghost_recv_[offr]),xgsize,MPI_DOUBLE,rank_xl_,
                 tagr,MPI_COMM_WORLD,&(xreqs_[1]));
#else
       // recv from left neighbor
       std::copy(xghost_send_+offr,xghost_send_+offr+xgsize,xghost_recv_+offl);
#endif
    }
    else{// left boundary is physical
       if(debug) fprintf(stderr,"rank=%d: xl is physical\n",rank_);
       // load boundary conditions to xghost_recv_
    }

    // right boundary
    if(rank_<rank_xr_ || strcmp(bound[1],"periodic")==0){// right boundary is MPI
       if(debug) fprintf(stderr,"rank=%d: xr is MPI\n",rank_);
#if USE_MPI
       // send to right neighbor
       MPI_Isend(&(xghost_send_[offr]),xgsize,MPI_DOUBLE,rank_xr_,
                 tagr,MPI_COMM_WORLD,&(xreqs_[2]));
       // recv from right neighbor
       MPI_Irecv(&(xghost_recv_[offr]),xgsize,MPI_DOUBLE,rank_xr_,
                 tagl,MPI_COMM_WORLD,&(xreqs_[3]));
#else
       // recv from right neighbor
       std::copy(xghost_send_+offl,xghost_send_+offl+xgsize,xghost_recv_+offr);
#endif
    }
    else{// left boundary is physical
       if(debug) fprintf(stderr,"rank=%d: xr is physical\n",rank_);
       // load boundary conditions to xghost_recv_
    }
 
    // wait for MPI communications
#if USE_MPI
    MPI_Waitall(4,xreqs_,xstats_);
#endif

    // load ghost cells
    if(debug>2){
       fprintf(stderr,"rank=%d:checking revc\n",rank_);
       checkMPI("xrecv.dat",xghost_send_,2*xgsize);
       fprintf(stderr,"rank=%d:finished checking recv\n",rank_);
    }
    grids->setGhostVec(-1,&(xghost_recv_[offl]),sendID); // left
    grids->setGhostVec(1,&(xghost_recv_[offr]),sendID); // right
    if(debug) fprintf(stderr,"rank=%d:Finished loading x Ghosts!\n",rank_);

    /* y field boundaries *********************************/
    // yz field boundaries 
    grids->updatePeriodicGhostCells();

}
