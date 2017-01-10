#if USE_MPI

#include <iostream>
#include <cmath>
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
        int dest_; // MPI destination send to, and recv from
        int stag_,rtag_; // MPI tags for send and recv
        double *send_;
        double *recv_;  
        Grid *grids_;     
};

BC_F_MPI::BC_F_MPI(int side, Domain* domain, Grid *grids, Input_Info_t *info){
 
    // load side 
    side_ = side; // inherited from BC_Field
    assert((side_>0 && side_<=3)||(side_<0 && side_ >=-3));
    if(debug)cerr<<"rank="<<rank_MPI<<":boundary side "<<side<< " is MPI\n";

    // determin MPI rank of neighbour
    int index = sideToindex(side);
    int *neighbours = domain->getNeighbours();
    dest_ = neighbours[index];
    if(debug)cerr<<"rank="<<rank_MPI<<":boundary side "<<side<< " neighbour is "<<dest_<<endl;

    // encode MPI send tag and recv tag in format tag = (sender,receiver)
    stag_ = size_MPI*rank_MPI + dest_;
    rtag_ = size_MPI*dest_ + rank_MPI;
    if(debug)cerr<<"rank="<<rank_MPI
                 <<": boudary side "<<side
                 << " stag="<< stag_ 
                 << ", rtag="<< rtag_<< endl;

    // load grids 
    grids_ = grids;
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


int BC_F_MPI::completeBC(int fieldID, int option){


    int size = grids_->getGhostVecSize(fieldID); 

    send_ = new double[size]; 
    recv_ = new double[size]; 

    MPI_Request reqs[2];
    MPI_Status stats[2];
    int err; 
    char fname[20]; // for debug file 
 
    // load ghost value to send 
    grids_->getGhostVec(side_,send_,fieldID);     

    if(debug>2){
       fprintf(stderr,"rank=%d:checking field send\n",rank_MPI);
       sprintf(fname,"field%d_send%d.dat",fieldID,side_);
       checkMPI(fname,send_,size);
       fprintf(stderr,"rank=%d:finished checking field send\n",rank_MPI);
    }
 
    // non-blocking send to neighbour
    err = MPI_Isend(send_,size,MPI_DOUBLE,dest_,stag_,MPI_COMM_WORLD,&reqs[0]);
    if(err) fprintf(stderr, "rank=%d: MPI_Isend field on error %d\n",rank_MPI,err);

    // non-blocking blocking recv from neighbor
    err = MPI_Irecv(recv_,size,MPI_DOUBLE,dest_,rtag_,MPI_COMM_WORLD,&reqs[1]);
    if(err) fprintf(stderr, "rank=%d: MPI_IRecv field on error %d\n",rank_MPI,err);

    // wait for communication
    err = MPI_Waitall(2,reqs,stats);
    if(err)fprintf(stderr, "rank=%d (1) MPI_Waitall error = %d\n",rank_MPI,err);
    if(debug>2) fprintf(stderr,"rank=%d: field sent and received on side %d\n",rank_MPI,side_);

    if(debug>2){
       fprintf(stderr,"rank=%d:checking field revc\n",rank_MPI);
       sprintf(fname,"field%d_recv%d.dat",fieldID,side_);
       checkMPI(fname,recv_,size);
       fprintf(stderr,"rank=%d:finished checking field recv\n",rank_MPI);
    }

    // unload received value to ghost 
    grids_->setGhostVec(side_,recv_,fieldID,option);   
 
    delete [] send_;
    delete [] recv_;

    return 0;
}

// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("MPI", makeBCField<BC_F_MPI>);
#endif
