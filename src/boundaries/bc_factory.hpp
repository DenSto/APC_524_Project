#ifndef BC_FACTORY
#define BC_FACTORY

#include <stdio.h>
#include <map>
#include <vector>
#include <sstream>
#include <stdexcept>
#include "boundary_particles.hpp"
#include "../domain/domain.hpp"
#include "assert.h"





/*!
 * A singleton class to handle registration of particle boundaries
 */
class BC_Factory {
public:
    typedef BC_Particle *(*Factory)(Domain* domain, int dim_Index, short isLeft, std::string type);
  /*  
   * Return the (singleton) BC_Factory object
   */
  static BC_Factory& getInstance() {
    static BC_Factory instance; // Guaranteed to be destroyed; Instantiated on first use.
    return instance;
  }

   /*
   	* Construct the boundary condition array (must be freed!)
	* Takes in an array of size 6.
    */
     //BC_Particle** constructConditions(Domain* domain, const std::string* types){
     BC_Particle** constructConditions(Domain* domain, const char (*bound)[32]){
                // convert c string for MPI to std::string for BC_Particle
                // size 32 is in correspondence definition in Input_Info_t
		std::string types[6];
		for(int i=0;i<6;i++){
		    types[i].assign(bound[i]);
		    //std::cerr << types[i] <<std::endl;
		}
#if USE_MPI
		std::string periodic ("periodic");
		std::string mpi ("MPI");
		int* nProc = domain->getnProcxyz();	
		int* myLoc = domain->getmyijk();	
#endif
		BC_Particle** ret = new BC_Particle*[6];
		assert(ret != NULL);
		for(int i = 0; i < 3; i++){
#if USE_MPI
			if(nProc[i] > 1){ // needs MPI
				int partitionIndex = myLoc[i];
				short inMiddle = (partitionIndex != 0) && (partitionIndex != nProc[i] - 1);
				// compare return 0 when equal
				if(periodic.compare(types[2*i])==0 || inMiddle){ // Two MPI communication loops
					// Left condition
					ret[2*i]=lookup(mpi)(domain,i,1,mpi);

					// Right condition
					ret[2*i+1]=lookup(mpi)(domain,i,0,mpi);
				} else { // One MPI communication loop
					// Physical side (to be calculated first)
					int MPIisLeft = 0;	
					if(partitionIndex == 0){ // Physical boundary on left
						ret[2*i]=lookup(types[2*i])(domain,i,1,types[2*i]);
						MPIisLeft=0;
					} else { // Physical boundary on right
						ret[2*i+1]=lookup(types[2*i+1])(domain,i,0,types[2*i+1]);
						MPIisLeft=1;
					}
					// MPI side (left or right, compute second)
					ret[2*i+1-MPIisLeft]=lookup(mpi)(domain,i,MPIisLeft, mpi);
				}
 			} else // treat serially
#endif
			{
				// Left condition
				ret[2*i]=lookup(types[2*i])(domain,i,1,types[2*i]);

				// Right condition
				ret[2*i+1]=lookup(types[2*i+1])(domain,i,0,types[2*i+1]);
			}
		}			 
		return ret;
	}

  /*  
   * Declare an particle boundary by type
   */
  void declare(const std::string &type, // the type of particle
               Factory factory             // a factory function to make the desired integrator
              ) {
    //fprintf(stderr,"call declare\n");
    registry_.insert(std::make_pair(type, factory));
  }

  /*
   * return the named factory function
   */
  Factory lookup(const std::string &type // the type in question
                ) {
    MapType::iterator i = registry_.find(type);
    if (i == registry_.end()) {     // not found
      std::ostringstream os;
      os << "Unable to find particle boundary " << type;
      throw std::runtime_error(os.str());
    }     

    return i->second;
  }

  /*
   * Return the sorted list of types of known boundary conditions
   */
  std::vector<const std::string *> types() const {
    std::vector<const std::string *> types;
    types.reserve(registry_.size());

    for (MapType::const_iterator ptr = registry_.begin(); ptr != registry_.end(); ++ptr) {
      types.push_back(&ptr->first);
    }     

    return types;
  }
private:
  BC_Factory() {}
  ~BC_Factory() {} 
  // Not implemented --- can't be called
  BC_Factory(const BC_Factory&);
  void operator=(const BC_Factory&);

  typedef std::map<std::string, Factory> MapType;
  MapType registry_; // BC_particle registry  registry
};


/*
 * An object which, when instantiated, registers a particle boundary condition
 */
struct RegisterParticleBoundary {
  RegisterParticleBoundary(const std::string &type, // type of boundary
                 BC_Factory::Factory factory   // factory function for integrator
                )
    {
      BC_Factory::getInstance().declare(type, factory);
    }
};

/*
 * A factory function for particle boudaries
 */
template<typename T>
BC_Particle *makeBCParticle(Domain* domain, int dim_Index, short isLeft, std::string type) {
    return new T(domain,dim_Index,isLeft, type);
}



#endif
