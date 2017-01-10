#ifndef PART_BC_FACTORY
#define PART_BC_FACTORY

#include <stdio.h>
#include <map>
#include <vector>
#include <sstream>
#include <stdexcept>
#include "particles_boundary.hpp"
#include "../IO/input.hpp"
#include "../domain/domain.hpp"
#include "assert.h"





/*!
 * A singleton class to handle registration of particle boundaries
 */
class Part_BC_Factory {
public:
    typedef BC_Particle *(*Factory)(Domain* domain, int dim_Index, short isRight, std::string type);
  /*  
   * Return the (singleton) Part_BC_Factory object
   */
  static Part_BC_Factory& getInstance() {
    static Part_BC_Factory instance; // Guaranteed to be destroyed; Instantiated on first use.
    return instance;
  }


  BC_Particle** constructConditions(Domain* domain, Input_Info_t* info);

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


  Input_Info_t* getInfo(){return info_;};
  void setInfo(Input_Info_t* info){info_ = info;};
  
private:
  Part_BC_Factory() {}
  ~Part_BC_Factory() {} 
  // Not implemented --- can't be called
  Part_BC_Factory(const Part_BC_Factory&);
  void operator=(const Part_BC_Factory&);

  typedef std::map<std::string, Factory> MapType;
  MapType registry_;  // BC_particle registry  registry
  Input_Info_t* info_; // BC_particle registry  info on particle tyes
};


/*
 * An object which, when instantiated, registers a particle boundary condition
 */
struct RegisterParticleBoundary {
  RegisterParticleBoundary(const std::string &type, // type of boundary
                 Part_BC_Factory::Factory factory   // factory function for integrator
                )
    {
      Part_BC_Factory::getInstance().declare(type, factory);
    }
};

/*
 * A factory function for particle boudaries
 */
template<typename T>
BC_Particle *makeBCParticle(Domain* domain, int dim_Index, short isRight, std::string type) {
    return new T(domain,dim_Index,isRight, type);
}



#endif
