#if USE_MPI

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string.h>
#include "mpi.h"
#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"
#include "../../IO/output.hpp"

using namespace std;

class BC_F_MPI : public BC_Field {
    public:
	BC_F_MPI(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_MPI();
        int sideToindex(int side);
	int completeBC(int sendID, int option);
    private:
        short same_; //same=+1: send and recv from the same side
                     //same=-1: send and recv from opposite sides
        int sendRank_; // MPI destination send to
        int recvRank_; // MPI destination recv from
        int stag_,rtag_; // MPI tags for send and recv
        double *sendBuf_;// send buffer
        double *recvBuf_;// recv buffer  
        double *ghostBuf_; // ghost buffer
        Grid *grids_;     
};

BC_F_MPI::BC_F_MPI(int side, Domain* domain, Grid *grids, Input_Info_t *info){
 
    // load grids 
    grids_ = grids;

    // load side 
    side_ = side; // inherited from BC_Field
    int dim = (int)(abs((double)side_)-1); 
    assert(dim>=0 && dim<3);
    //if(debug)cerr<<"rank="<<rank_MPI<<": boundary side "<<side_<< " is MPI\n";
    if(debug)fprintf(stderr,"rank=%d: boundary side %d is MPI\n",rank_MPI,side_);

    // if in the middle or periodic, send this side, recv the other side
    // otherwise, send and recv from the same MPI side
    int ind_this = sideToindex(side);
    int ind_other= sideToindex(-side);
    bool isPeriodic = (strcmp(info->fields_bound[ind_this],"periodic")==0);
    //if(debug>1&&isPeriodic)cerr<<"rank="<<rank_MPI<<": boundary side is Periodic\n";
    if(debug>1&&isPeriodic)fprintf(stderr,"rank=%d: side %d is Periodic\n",rank_MPI,side_);

    // determine MPI locations and neighbours 
    int* nProc = domain->getnProcxyz();
    int* myLoc = domain->getmyijk();
    int* neigh = domain->getNeighbours();

    int partitionIndex = myLoc[dim];
    short inMiddle = (partitionIndex != 0 && (partitionIndex != nProc[dim]- 1));

    int rank_this = neigh[ind_this]; // MPI neighbour on this side
    int rank_other= neigh[ind_other];// MPI neighbour on the other side    

    // determine which neighbour to send and recv
    // handled both to avoid MPI deadlock
    if(inMiddle || isPeriodic){
        sendRank_ = rank_this; //send to neighbour on this side
        recvRank_ = rank_other;//recv from neighbour on the other side
        same_=-1;
    } else {
        sendRank_ = rank_this;//send to this side
        recvRank_ = rank_this;//recv from this side
        same_=1;
    }
    //if(debug)cerr<<"rank="<<rank_MPI<<": boundary side "<<side_
    //             << " sendRank="<<sendRank_
    //             << ", recvRank="<<recvRank_<<endl;
    if(debug)fprintf(stderr,"rank=%d: boundary side=%d, sendRank=%d, recvRank=%d\n",
                             rank_MPI,side,sendRank_,recvRank_);


    // encode MPI send tag and recv tag in format tag = (sender,receiver)
    stag_ = size_MPI*rank_MPI + sendRank_;
    rtag_ = size_MPI*recvRank_+ rank_MPI;
    //if(debug)cerr<<"rank="<<rank_MPI
    //             <<": boudary side "<<side_
    //             << " stag="<< stag_ 
    //             << ", rtag="<< rtag_<< endl;
    if(debug)fprintf(stderr,"rank=%d: boudary side=%d, stag=%d, rtag=%d\n",
                            rank_MPI,side_,stag_,rtag_);
}

BC_F_MPI::~BC_F_MPI(){
}

//! Function convert side (-3,-2,-1,1,2,3) to index (4,2,0,1,3,5) 
int BC_F_MPI::sideToindex(int side){
    double index;
    index = 2*abs((double)side+0.25)-1.5;
    return (int)index;
}
/* test sideToindex
int sides[6] = {-3,-2,-1,1,2,3};
for(int i=0;i<6;i++){
    cerr<<"side="<<sides[i]<<" -> index="<<sideToindex(sides[i])<< endl;
}*/


//! complete MPI field boundary condition
/*! option = 0: load physical on this side, and send to neighbour 
 *              recv from neighbor, use it to replace ghost on this side
 *  option = 1: take ghost on this side, and add to physical on this side
 *              no MPI communication is needed, only local operations
 */             
int BC_F_MPI::completeBC(int fieldID, int option){

    assert(option==0 || option==1);

    int size = grids_->getGhostVecSize(fieldID); 

    if(option==0){//send physical, recv to ghost
        if(debug>1)fprintf(stderr,"rank=%d: side=%d MPI for field %d...\n",
                                   rank_MPI,side_,fieldID);

        sendBuf_ = new double[size]; 
        recvBuf_ = new double[size]; 

        MPI_Request req;
        int err; 
        char fname[20]; // for debug file 
    
        // load ghost value to send. Always send this side of ghost
        grids_->getGhostVec(side_,sendBuf_,fieldID,0);//0: get physical
   
        if(debug>2){
           fprintf(stderr,"rank=%d:checking field send\n",rank_MPI);
           sprintf(fname,"field%d_send%d.dat",fieldID,side_);
           checkMPI(fname,sendBuf_,size);
           fprintf(stderr,"rank=%d:finished checking field send\n",rank_MPI);
        }
    
        // non-blocking send to neighbour
        err=MPI_Isend(sendBuf_,size,MPI_DOUBLE,sendRank_,stag_,MPI_COMM_WORLD,&req);
        if(err) fprintf(stderr, "rank=%d: MPI_Isend field on error %d\n",rank_MPI,err);
    
        // blocking recv from neighbor, the nieghbour should have already sent 
        err=MPI_Recv(recvBuf_,size,MPI_DOUBLE,recvRank_,rtag_,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if(err) fprintf(stderr, "rank=%d: MPI_Recv field on error %d\n",rank_MPI,err);
    
        // wait for communication
        err=MPI_Wait(&req,MPI_STATUS_IGNORE);
        if(err) fprintf(stderr, "rank=%d MPI_Wait error = %d\n",rank_MPI,err);
        if(debug>2) fprintf(stderr,"rank=%d: field received on side %d\n",rank_MPI,side_);
    
        if(debug>2){
           fprintf(stderr,"rank=%d:checking field revc\n",rank_MPI);
           sprintf(fname,"field%d_recv%d.dat",fieldID,side_);
           checkMPI(fname,recvBuf_,size);
           fprintf(stderr,"rank=%d:finished checking field recv\n",rank_MPI);
        }
    
        // unload received value to ghost, unload to this or oppisite side
        // same side: same =+1. Opposite side: same =-1
        grids_->setGhostVec(same_*side_,recvBuf_,fieldID,0);//0: replace ghost

        delete [] sendBuf_;
        delete [] recvBuf_;

    } else { // load ghost, sum to physical, this side only 
 
        if(debug>1)fprintf(stderr,"rank=%d: side=%d MPI local sum for field %d...\n",
                                   rank_MPI,side_,fieldID);
        ghostBuf_ = new double[size]; 
        // load ghost on this side
        grids_->getGhostVec(side_,ghostBuf_,fieldID,1); //1: get ghost
        // sum to physical on this side
        grids_->setGhostVec(side_,ghostBuf_,fieldID,1); //1: sum to physical
        delete [] ghostBuf_;

    }

    return 0;
}

// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("MPI", makeBCField<BC_F_MPI>);
#endif
