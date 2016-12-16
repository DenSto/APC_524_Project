#include <stdio.h>
#include <assert.h>
#include "domain.hpp"
#include "../IO/output.hpp"

#if USE_MPI
#include "mpi.h"
#endif

//! Pass fields across MPI boundaries, or execute physical boundary conditions
void Domain::PassFields(Grid *grids, Input_Info_t *input_info){
    int tagl=0; // MPI tag: to left or from right
    int tagr=1; // MPI tag: to right or from left
    char (*bound)[32] = input_info->fields_bound;
    
    /* x field boundaries *********************************/
    int xgsize = grids->getGhostVecSize(); // ghost size in x direction
    assert(xgsize>0);   
    fprintf(stderr,"rank=%d: xgsize=%d\n",rank_,xgsize);

    // load buffer xghost_send_
    grids->getGhostVec(-1, xghost_send_+tagl*xgsize);// left
    grids->getGhostVec(1, xghost_send_+tagr*xgsize);// right
    fprintf(stderr,"rank=%d:checking send\n",rank_);
    checkMPI("./xsend.dat",xghost_send_,2*xgsize);
    fprintf(stderr,"rank=%d:finished checking send\n",rank_);

    // left boundary
    if(rank_>rank_xl_ || strcmp(bound[0],"periodic")==0){// left boundary is MPI
       fprintf(stderr,"rank=%d: xl is MPI\n",rank_);
#if USE_MPI
       // send to left neighbor
       // shift pointer of buffer
       MPI_Isend(xghost_send_+tagl*xgsize,xgsize,MPI_DOUBLE,rank_xl_,
                 tagl,MPI_COMM_WORLD,&(xreqs_[0]));
       // recv from left neighbor
       // shift pointer of buffer
       MPI_Irecv(xghost_recv_+tagl*xgsize,xgsize,MPI_DOUBLE,rank_xl_,
                 tagr,MPI_COMM_WORLD,&(xreqs_[1]));
#else
       // recv from left neighbor
       //std::copy(xghost_send_+xgsize,xghost_send_+2*xgsize,xghost_recv_);
#endif
    }
    else{// left boundary is physical
       fprintf(stderr,"rank=%d: xl is physical\n",rank_);
       // load boundary conditions to xghost_recv_
    }

    // right boundary
    if(rank_<rank_xr_ || strcmp(bound[1],"periodic")==0){// right boundary is MPI
       //fprintf(stderr,"rank=%d: xr is MPI\n",rank_);
#if USE_MPI
       // send to right neighbor
       // shift pointer of buffer
       MPI_Isend(xghost_send_+tagr*xgsize,xgsize,MPI_DOUBLE,rank_xr_,
                 tagr,MPI_COMM_WORLD,&(xreqs_[2]));
       // recv from right neighbor
       // shift pointer of buffer
       MPI_Irecv(xghost_recv_+tagr*xgsize,xgsize,MPI_DOUBLE,rank_xr_,
                 tagl,MPI_COMM_WORLD,&(xreqs_[3]));
#else
       // recv from right neighbor
       //std::copy(xghost_send_,xghost_send_+xgsize,xghost_recv_+xgsize);
#endif
    }
    else{// left boundary is physical
       //fprintf(stderr,"rank=%d: xr is physical\n",rank_);
       // load boundary conditions to xghost_recv_
    }
 
    // wait for MPI communications
#if USE_MPI
    MPI_Waitall(4,xreqs_,xstats_);
#endif

    // load ghost cells
    checkMPI("./xrecv.dat",xghost_send_,2*xgsize);
    grids->setGhostVec(-1,&xghost_recv_[tagl*xgsize]); // left
    grids->setGhostVec(1,&xghost_recv_[tagr*xgsize]); // right
    fprintf(stderr,"rank=%d:Finished loading x Ghosts!\n",rank_);

    // yz field boundaries 
    grids->updatePeriodicGhostCells();

}
