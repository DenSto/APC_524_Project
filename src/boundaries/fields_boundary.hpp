#ifndef BOUNDARY_FIELDS_HPP
#define BOUNDARY_FIELDS_HPP

#include "../grid/grid.hpp"
#include <string>

//! Class which defines a field boundary condition
class BC_Field {
	public:
		virtual ~BC_Field() {};
	private:
		virtual int completeBC(void) = 0; // Complete boundary conditions.
		std::string type_;
};

#endif
