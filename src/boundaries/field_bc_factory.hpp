#ifndef FIELD_BC_FACTORY
#define FIELD_BC_FACTORY

#include <stdio.h>
#include <map>
#include <vector>
#include <sstream>
#include <stdexcept>
#include "fields_boundary.hpp"
#include "../IO/input.hpp"
#include "../domain/domain.hpp"
#include "../grid/grid.hpp"
#include "assert.h"

//! A singleton class to handle registration of field boundaries/
class Field_BC_Factory {
    public:
        // Return the (singleton) Field_BC_Factory object
        static Field_BC_Factory& getInstance() {
            static Field_BC_Factory instance; 
            return instance;
        }
       
        void constructConditions(Domain* domain, Grid *grids, const char (*bound)[NCHAR]);
       
        typedef BC_Field *(*Factory)(Domain* domain, Grid *grids, int side);
       
        // Declare an particle boundary by type
        void declare(const std::string &type, // the type of particle
                     Factory factory          // a factory function to make the desired boundary
                    ) {
          //fprintf(stderr,"call declare\n");
          registry_.insert(std::make_pair(type, factory));
        }
       
        // return the named factory function
        Factory lookup(const std::string &type) {
            MapType::iterator i = registry_.find(type);
            if (i == registry_.end()) {     // not found
              std::ostringstream os;
              os << "Unable to find field boundary " << type;
              throw std::runtime_error(os.str());
            }     
            return i->second;
        }
       
        //Return the sorted list of types of known boundary conditions
        std::vector<const std::string *> types() const {
            std::vector<const std::string *> types;
            types.reserve(registry_.size());
       
            for(MapType::const_iterator ptr=registry_.begin();ptr!=registry_.end();++ptr){
              types.push_back(&ptr->first);
            }     
            return types;
        }

    private:
        Field_BC_Factory() {}
        ~Field_BC_Factory() {} 
        void operator=(const Field_BC_Factory&);
       
        typedef std::map<std::string, Factory> MapType;
        MapType registry_; // BC_field registry

};


//! An object which, when instantiated, registers a field boundary condition
struct RegisterFieldBoundary {
  RegisterFieldBoundary(const std::string &type, // type of boundary
                 Field_BC_Factory::Factory factory   // factory function for integrator
                )
    {
      Field_BC_Factory::getInstance().declare(type, factory);
    }
};

//! A factory function for particle boudaries
template<typename T>
BC_Field *makeBCField(Domain* domain, Grid *grids, int side) {
    return new T(domain,grids,side);
}

#endif
