#include <string>
#include <math.h>
#include "fieldBC.hpp"
#include "grid.hpp"


FieldBC::FieldBC( std::string &fieldStr, int dim, bool side, double amp, double omega, double phase):
   fieldStr_(fieldStr)
{
   dim_ = dim;
   side_ = side;
   amp_ = amp;
   omega_ = omega;
   phase_ = phase;
};

//! Apply boundary condition to grid
/*!
   Uses setFieldAlongEdge method in grid to add field to grid.
*/
void FieldBC::applyBCs( double t, double dt, Grid &grid) {
   double fieldVal = amp_ * cos( omega_ * t + phase_ );
   grid.setFieldAlongEdge( fieldStr_, dim_, side_, fieldVal);
};
