#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_P_MPI : public BC_Particle {
	public:
		BC_P_MPI(Domain* domain, int dim_Index, std::string type);
		~BC_P_MPI();
		void computeParticleBCs(std::vector<Particle> pl);
		void completeBC();
	private:
		int particle_BC(double *x, double *v, double xMin, double xMax);
		double xMin_;
		double xMax_;
		int dim_index_;
		std::string type_;

		short isLeft_;
		int targetRank_;
		double lengthShift_;
};


BC_P_MPI::BC_P_MPI(Domain* domain, int dim_Index, std::string type)
	:	dim_index_(dim_Index),
		type_(type)
		{}

BC_P_MPI::~BC_P_MPI(){
}

void BC_P_MPI::completeBC(){
	//Send/Receives go here!
}

int BC_P_MPI::particle_BC(double* x, double* v, double xMin, double xMax){
	// Load particles into buffer if needed
	return 0;
}
