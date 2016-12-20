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


  BC_Particle** constructConditions(Domain* domain, const char (*bound)[32]);

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
