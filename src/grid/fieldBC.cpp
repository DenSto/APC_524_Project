
#include <string>
#include <math.h>
#include "fieldBC.hpp"
#include "grid.hpp"


FieldBC::FieldBC( std::string &fieldStr, int dim, bool edge, double amp, double omega, double phase):
   fieldStr_(fieldStr),
   dim_(dim),
   edge_(edge),
   amp_(amp),
   omega_(omega),
   phase_(phase)
{};

//! Apply boundary condition to grid
/*!
   Uses setFieldAlongEdge method in grid to add field to 
*/
void FieldBC::applyBCs( double t, Grid &grid) {
   double fieldVal = amp_ * cos( omega_ * t + phase_ );
   grid.setFieldAlongEdge( fieldStr_, dim_, edge_, fieldVal);
};