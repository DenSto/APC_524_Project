#ifndef BOUNDARY_FIELDS_HPP
#define BOUNDARY_FIELDS_HPP

#include <string>

//! Class which defines a field boundary condition
class BC_Field {
    public:
	virtual ~BC_Field() {};
	virtual int completeBC(int sendID, int option) = 0; // Complete boundary conditions.
        // sendID determine which field to complete, see ../grid/spookyGrid.cpp for details
        // option = 0: replace, 1: sum

    protected:
        int side_;
};

#endif
