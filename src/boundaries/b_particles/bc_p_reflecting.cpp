#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_P_Reflecting : public BC_Particle {
	public:
		BC_P_Reflecting(Domain* domain, int dim_Index, std::string type);
		~BC_P_Reflecting();
		void computeParticleBCs(std::vector<Particle> pl);
		void completeBC();
	private:
		int particle_BC(double *x, double *v, double xMin, double xMax);
		double xMin_;
		double xMax_;
		int dim_index_;
		std::string type_;
};

BC_P_Reflecting::BC_P_Reflecting(Domain* domain, int dim_Index, std::string type) 
	:	dim_index_(dim_Index),
		type_(type)
		{}

BC_P_Reflecting::~BC_P_Reflecting(){
}

void BC_P_Reflecting::completeBC(){}

int BC_P_Reflecting::particle_BC(double* x, double* v, double xMin, double xMax){
	if(*x > xMax){
		*x = 2.0*xMax - *x;
		*v=-*v;
	} 
	if(*x < xMin){
		*x = 2.0*xMin - *x;
		*v=-*v;
	}
	return 0;
}

static RegisterParticleBoundary instance("reflecting", makeBCParticle<BC_P_Reflecting>);
