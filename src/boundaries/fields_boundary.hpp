#ifndef BOUNDARY_FIELDS_HPP
#define BOUNDARY_FIELDS_HPP

#include <string>

//! Class which defines a field boundary condition
class BC_Field {
    public:
	virtual ~BC_Field() {};
	virtual int completeBC(void) = 0; // Complete boundary conditions.

    protected:
        int side_;

};

#endif
