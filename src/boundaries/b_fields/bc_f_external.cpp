#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_F_External : public BC_Field {
	public:
		BC_F_External(Domain* domain, Grid *grids, int side);
		~BC_F_External();
		int completeBC();
//		void computeParticleBCs(std::vector<Particle> *pl);
	private:
//		int particle_BC(Particle* p);
};

BC_F_External::~BC_F_External(){
}

BC_F_External::BC_F_External(Domain* domain, Grid *gids, int side){
/*
	xMin_ = domain->getxyz0()[dim_index_];
        xMax_ = xMin_+domain->getLxyz()[dim_index_];
	if(debug>1)fprintf(stderr,"rank=%d:dim=%d,isRight=%d,reflect_BC,xMin=%f,xMax=%f\n",
                                   rank_MPI,dim_index_,isRight_,xMin_,xMax_); 	

*/
}

int BC_F_External::completeBC(){
	// Nothing to do
	return 0;
}

/*
int BC_P_Reflecting::particle_BC(Particle* p){
	if(p->x[dim_index_] > xMax_ && isRight_){
		p->x[dim_index_] = 2.0*xMax_ - p->x[dim_index_];
		p->v[dim_index_]=-p->v[dim_index_];
	} 
	if(p->x[dim_index_] < xMin_ && !isRight_){
		p->x[dim_index_] = 2.0*xMin_ - p->x[dim_index_];
		p->v[dim_index_]=-p->v[dim_index_];
	}
	return 0;
}
*/

// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("external", makeBCField<BC_F_External>);
