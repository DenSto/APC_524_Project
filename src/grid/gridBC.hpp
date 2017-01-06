#ifndef GRID_BC_HPP
#define GRID_BC_HPP

#include "grid.hpp"

//! Abstract class for supplying boundary conditions to field grid.
/*!
   side = -1: x left, side = +1: x right \n
   side = -2: y left, side = +2: y right \n
   side = -3: z left, side = +3: y right \n

   dim_ = abs(side_)-1
   dim_ = 0: x boundaries
   dim_ = 1: y boundaries
   dim_ = 2: z boundaries
*/
class GridBC {
public:
   virtual ~GridBC() {};
   // t: current time, dt: for phase shift
   virtual void applyBCs( double t, double dt, Grid *grids) = 0;

protected:
   int side_;
   int dim_; 
};


#endif
