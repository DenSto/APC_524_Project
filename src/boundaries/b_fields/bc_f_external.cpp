#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"

class BC_F_External : public BC_Field {
    public:
	BC_F_External(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_External();
	int completeBC();
    private:
        Grid *grids_;
        Input_Info_t *input_info_;
};

BC_F_External::BC_F_External(int side, Domain* domain, Grid *grids, Input_Info_t *info){
    side_ = side;
    assert((side_>0 && side_<=3)||(side_<0 && side_ >=-3));
    if(debug)fprintf(stderr,"rank=%d:boundary side %d is external\n",rank_MPI,side);

    grids_ = grids;
    
    input_info_ = info;
}

BC_F_External::~BC_F_External(){
}

int BC_F_External::completeBC(){
	// Nothing to do
	return 0;
}

// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("external", makeBCField<BC_F_External>);
