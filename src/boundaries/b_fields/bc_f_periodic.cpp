#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"

class BC_F_Periodic : public BC_Field {
    public:
        BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_Periodic();
	int completeBC();
    private:
        Grid *grids_;     
};

BC_F_Periodic::BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info){
    side_ = side;
    assert((side_>0 && side_<=3)||(side_<0 && side_ >=-3));
    if(debug)fprintf(stderr,"rank=%d:boundary side %d is periodic\n",rank_MPI,side);

    grids_ = grids;
    ghostTmp_ = grids_->ghostTmp_; 
}


BC_F_Periodic::~BC_F_Periodic(){
}


int BC_F_Periodic::completeBC(){

    // load the opposite side to tmp
    grids_->getGhostVec(-side_,ghostTmp_,-1);
    // unload this side fom tmp
    grids_->setGhostVec(side_,ghostTmp_,-1);

    return 0;
}


// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("periodic", makeBCField<BC_F_Periodic>);
