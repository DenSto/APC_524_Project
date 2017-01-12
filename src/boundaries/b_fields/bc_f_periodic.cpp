#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"

class BC_F_Periodic : public BC_Field {
    public:
        BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_Periodic();
	int completeBC(int fieldID, int option);
    private:
        double* ghost_;
        Grid *grids_;     
};

BC_F_Periodic::BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info){

    side_ = side; // inherited from BC_Field
    assert((side_>0 && side_<=3)||(side_<0 && side_ >=-3));
    if(debug)fprintf(stderr,"rank=%d:boundary side %d is periodic\n",rank_MPI,side);

    grids_ = grids;

}


BC_F_Periodic::~BC_F_Periodic(){
}

//! complete periodic field boundary condition
/*! option = 0: load physical of the other side
 *              replace ghost on this side
 *  option = 1: load ghost on this side
 *              sum ghost to physical on this side
 */
int BC_F_Periodic::completeBC(int fieldID, int option){

    int size;
    size = grids_->getGhostVecSize(fieldID); 
    ghost_ = new double[size]; 

    assert(option==0 || option==1);
    if(option==0){
        // load physical on the other side 
        grids_->getGhostVec(-side_,ghost_,fieldID,0);
        // replace ghost on this side
        grids_->setGhostVec(side_,ghost_,fieldID,0);
    } else {
        // load ghost on this side  
        grids_->getGhostVec(side_,ghost_,fieldID,1);
        // sum ghost to physical on this side
        grids_->setGhostVec(side_,ghost_,fieldID,1);
    }

    delete [] ghost_;

    return 0;
}


// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("periodic", makeBCField<BC_F_Periodic>);
