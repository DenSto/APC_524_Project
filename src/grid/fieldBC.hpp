#ifndef FIELD_BC_HPP
#define FIELD_BC_HPP

#include "gridBC.hpp"
#include <string>

class Grid;


//! Class for supplying boundary conditions in a single field to field grid.
/*!
   Boundary conditions are of form: \n
   amp * cos( omega * t + phase) \n
   along plane perpendicular to dimension dim (0 = x, 1 = y, 2 = z)
   on edge ( false = left, true = right) \n
   fieldStr one of Ex, Ey, Ez, Bx, By, Bz
*/
class FieldBC : public GridBC {
public:
   FieldBC( std::string &fieldStr, int dim, bool edge, double amp, double omega, double phase);
   void applyBCs( double t, double dt, Grid &grid);

private:
   std::string fieldStr_;
};


#endif
