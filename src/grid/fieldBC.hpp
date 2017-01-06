#ifndef FIELD_BC_HPP
#define FIELD_BC_HPP

#include "gridBC.hpp"
#include "../IO/input.hpp"

//! Class for supplying boundary conditions in a single field to field grid.
/*!
   Boundary conditions are of linear superpositions of the form: \n
   peakamp * cos( omega * t + phase) * exp(-(t-t0)^2/w^2) \n
   along plane perpendicular to side 
   side = -1: x left, side = +1: x right \n
   side = -2: y left, side = +2: y right \n
   side = -3: z left, side = +3: y right \n

   fieldStr one of Ex, Ey, Ez, Bx, By, Bz
*/
class FieldBC : public GridBC {
public:
   FieldBC(int side, Input_Info_t *input_info);
   void applyBCs(double t, Grid *grids);

private:
   Input_Info_t *input_info_;
};


#endif
