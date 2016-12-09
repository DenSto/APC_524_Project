#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_P_Periodic : public BC_Particle {
	public:
		BC_P_Periodic(Domain* domain, int dim_Index, std::string type);
		~BC_P_Periodic();
		void computeParticleBCs(std::vector<Particle> pl);
		void completeBC();
	private:
		int particle_BC(double *x, double *v, double xMin, double xMax);
		double xMin_;
		double xMax_;
		int dim_index_;
		std::string type_;

		short isLeft_;
};

BC_P_Periodic::BC_P_Periodic(Domain* domain, int dim_Index, std::string type) 
	:	dim_index_(dim_Index),
		type_(type)
		{}

BC_P_Periodic::~BC_P_Periodic(){
}

void BC_P_Periodic::completeBC(){}

int BC_P_Periodic::particle_BC(double* x, double* v, double xMin, double xMax){
	if(*x > xMax)
		*x-= (xMax-xMin);
	if(*x < xMin)
		*x+= (xMax-xMin);
	return 0;
}

static RegisterParticleBoundary instance("periodic", makeBCParticle<BC_P_Periodic>);
