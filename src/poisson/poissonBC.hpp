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
class PoissonBC : public GridBC {
public:
   PoissonBC(int side, Input_Info_t *input_info);
   ~PoissonBC(void);
   void applyBCs(double fieldID, double option, Grid *grids);

private:
   int size_;
   int ifLoad_[4]; //=[ifphi,ifAx,ifAy,ifAz]
                  // ifphi=0: boundary condition have not been loaded
                  // ifphi=1: boundary condition have been loaded in u1
                  // ifphi>1: boundary condition have been loaded for
                  //          both u1 and u2, no need to load again.
   double *ghost_;
   double phiA_[4]; //=[phi,Ax,Ay,Az]
};


#endif
