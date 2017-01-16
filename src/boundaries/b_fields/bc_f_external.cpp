#include <string.h>
#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"
#include "../../grid/lightBC.hpp"
#include "../../poisson/poissonBC.hpp"

class BC_F_External : public BC_Field {
    public:
	BC_F_External(int side, Domain* domain, Grid *grids, Input_Info_t *info);
	~BC_F_External();
	int completeBC(int fieldID);
    private:
        Grid *grids_;
        LightBC *lightbc_; // inject transverse EM waves + background

        bool support_poisson_; // whether support poisson initialization
        PoissonBC *poissonbc_; // constant boundary for poisson initialization 
};

BC_F_External::BC_F_External(int side, Domain* domain, Grid *grids, Input_Info_t *info){
    side_ = side;
    assert((side_>0 && side_<=3)||(side_<0 && side_ >=-3));
    if(debug)fprintf(stderr,"rank=%d:boundary side %d is external\n",rank_MPI,side);

    grids_ = grids;

    // instantiate the class that handles the injection
    lightbc_ = new LightBC(side_,info);   

    // if initialize with poisson, than equipe poisson boundary conditions
    if(info->restart==0 && strcmp(info->fields_init,"poisson")==0){
        support_poisson_=true;
        poissonbc_ = new PoissonBC(side_,info,grids);
    }else{
        support_poisson_=false;
    }
}

BC_F_External::~BC_F_External(void){ 
   delete lightbc_;
   if(support_poisson_) delete poissonbc_;
}

int BC_F_External::completeBC(int fieldID){

    if(debug>1)fprintf(stderr,"    rank=%d: executing external boundary on side %d\n",
                    rank_MPI,side_); 
 
    if(fieldID==-1){//E,B fields 
        lightbc_->applyBCs(time_phys,dt_phys,grids_);	
    }else if(support_poisson_ && fieldID>0){
        //individual field for poisson solver
        poissonbc_->applyBCs(fieldID);
    }

    return 0;
}

// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("external", makeBCField<BC_F_External>);

