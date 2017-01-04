#ifndef LIGHT_BC_HPP
#define LIGHT_BC_HPP

#include "gridBC.hpp"


class Grid;


//! Class for supplying light-wave boundary conditions to field grid.
/*!
   Boundary conditions are of form: \n
   E = amp * cos( omega * t + phase) \n
   along plane perpendicular to dimension dim (0 = x, 1 = y, 2 = z)
   on edge ( false = left, true = right) \n
   fieldStr one of Ex, Ey, Ez, Bx, By, Bz
*/
class LightBC : public GridBC {
public:
   LightBC( int dim, bool edge, double amp, double omega, double phase);
   void applyBCs( double t, double dt, Grid &grid) ;

};


#endif