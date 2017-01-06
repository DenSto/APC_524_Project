#include <assert.h>
#include "globals.hpp"
#include "fieldBC.hpp"


FieldBC::FieldBC(int side, Input_Info_t *input_info):
{
   side_ = side;
   input_info_ = input_info;
   assert(input_info_!=NULL);
};

//! Apply boundary condition to grid
/*!
   Uses setFieldAlongEdge method in grid to add field to grid.
*/
void FieldBC::applyBCs(double t, Grid *grids) {

   // load input info relevant to boundary specified by side
   int nwaves_tot = input_info_->nwaves;
   int nwaves = 0;
   for(int i=0;i<nwaves_tot;i++){
      if(input_info_->inSide[i]==side){nwaves+=1;}
   }
   if(debug>1)fprintf(stderr,"rank=%d: side %d has %d waves injected\n",
                              rank_MPI,side_,nwaves);

   double fieldVal = amp_ * cos( omega_ * t + phase_ );
   grid.setFieldAlongEdge( fieldStr_, dim_, edge_, fieldVal);
};
