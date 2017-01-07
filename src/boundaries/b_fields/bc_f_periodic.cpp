#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"

class BC_F_Periodic : public BC_Field {
    public:
        BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_Periodic();
	int completeBC();
    private:
        double* ghostEB_;
        double *ghostRJ_;  
        Grid *grids_;     
};

BC_F_Periodic::BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info){

    side_ = side; // inherited from BC_Field
    assert((side_>0 && side_<=3)||(side_<0 && side_ >=-3));
    if(debug)fprintf(stderr,"rank=%d:boundary side %d is periodic\n",rank_MPI,side);

    grids_ = grids;

    int ghostsize;

    ghostsize = grids_->getGhostVecSize(-1); // E,B
    ghostEB_ = new double[ghostsize]; 

    ghostsize = grids_->getGhostVecSize(-2); // rho, J
    ghostRJ_ = new double[ghostsize];

}


BC_F_Periodic::~BC_F_Periodic(){
    delete [] ghostEB_;
    delete [] ghostRJ_;
}


int BC_F_Periodic::completeBC(){

    // load physical to tmp from the opposite side 
    grids_->getGhostVec(-side_,ghostEB_,-1);
    grids_->getGhostVec(-side_,ghostRJ_,-2);

    // unload tmp to Ghost on this side
    // Ghost does not overwrite physical on this side
    grids_->setGhostVec(side_,ghostEB_,-1,0); // replace(0) E,B(-1) in ghost
    grids_->setGhostVec(side_,ghostRJ_,-2,0); // replace(0) R,J(-2) in ghost 

    return 0;
}


// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("periodic", makeBCField<BC_F_Periodic>);
