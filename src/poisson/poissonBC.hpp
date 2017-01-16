#ifndef POISSON_BC_HPP
#define POISSON_BC_HPP

#include <vector>
#include <string>
#include "../grid/gridBC.hpp"
#include "../IO/input.hpp"
#include "poisson.hpp"


//! Class for supplying boundary conditions for poisson solver.
/*!
       side = -1: x left, side = +1: x right \n
       side = -2: y left, side = +2: y right \n
       side = -3: z left, side = +3: y right \n
*/
class PoissonBC {
public:
   PoissonBC(int side, Input_Info_t *input_info, Grid *grids);
   ~PoissonBC(void);
   void applyBCs(int fieldID);
   int fieldIDToIndex(int fieldID); 

private:
   Grid *grids_;

   double fieldVal_;
   double ***field_;
   double ****fieldPtr_;

   int dim_; //1:x, 2:y, 3:z
   int offset_; //where to place the boudary value
               // left:  first physical
               // right: last ghost
   int fid_;
   double phiA_[4]; 
};


#endif
