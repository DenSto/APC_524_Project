#ifndef GRID_BC_HPP
#define GRID_BC_HPP


#include <string>

class Grid;


//! Abstract lass for supplying boundary conditions to field grid.
/*!
   Boundary conditions are of form: \n
   amp * cos( omega * t + phase) \n
   along plane perpendicular to dimension dim (0 = x, 1 = y, 2 = z)
   on edge ( false = left, true = right) \n
*/
class GridBC {
public:
   virtual ~GridBC() {};
   virtual void applyBCs( double t, double dt, Grid &grid) = 0;

protected:
   int dim_;
   bool edge_;
   double amp_;
   double omega_;
   double phase_;

};


#endif