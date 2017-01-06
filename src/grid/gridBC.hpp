#ifndef GRID_BC_HPP
#define GRID_BC_HPP

#include "grid.hpp"

//! Abstract class for supplying boundary conditions to field grid.
/*!
   Boundary conditions are of form: \n
   peakamp * cos( omega * t + phase) * exp(-(t-delay)^2*invWidth^2) \n
   along plane perpendicular to side 
   side = -1: x left, side = +1: x right \n
   side = -2: y left, side = +2: y right \n
   side = -3: z left, side = +3: y right \n

*/
class GridBC {
public:
   virtual ~GridBC() {};
   virtual void applyBCs( double t, double dt, Grid *grids) = 0;

protected:
   int side_;
   int nwaves_;
   double *peakamps_;
   double *omegas_;
   double *phases_;
   double *invWidths_;
   double *delays_;

};


#endif
