#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"

class BC_F_Periodic : public BC_Field {
    public:
        BC_F_Periodic(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_Periodic();
	int completeBC(int fieldID);
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
/*!  
 *  Treat E,B,A,phi in one way
 *  Treat Rho, J in another way
 */
int BC_F_Periodic::completeBC(int fieldID){

    int size;
    size = grids_->getGhostVecSize(fieldID); 
    ghost_ = new double[size]; 

    if(fieldID==-2){ //Rho and orthogonal J's
        if(side_<0){// use left side to handle, right side do nothing 
            // load primary ghost on left 
            grids_->getGhostVec(-side_,ghost_,-2);
            // sum to physical on right
            grids_->setGhostVec(side_,ghost_,-2);
        }
    } else { // other fields
        // load physical on the other side  
        grids_->getGhostVec(-side_,ghost_,fieldID);
        // replace ghost on this side
        grids_->setGhostVec(side_,ghost_,fieldID);
    }

    delete [] ghost_;

    return 0;
}


// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("periodic", makeBCField<BC_F_Periodic>);
