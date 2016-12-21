#include <stdlib.h>
#include <stdio.h>
#include "interpolate.hpp"

Interpolator::Interpolator(){
}

Interpolator::~Interpolator(){
}

void Interpolator::interpolate_fields(double* pos, double* lcell, double* cellvars, Field_part* field) {
  int i,j,k;
  double ds[3]; /*projected distances from particles to "least" cell faces*/
  double eField, bField; /*field variables for calculation*/

  for (i=0; i < 3; i++) {
    ds[i] = pos[i] - cellvars[i];
  }

  for (i=0; i < 3; i++) {
    j = (i+1) % 3; /*Modular arithmetic cycles over Cartesian directions*/
    k = (i+2) % 3;

    eField = 0; /*Field additions are weighted by catercorner area.*/
    eField += (lcell[j]-ds[j]) * (lcell[k]-ds[k]) * cellvars[3+4*i]; /*1st edge*/
    eField += ds[j] * (lcell[k]-ds[k]) * cellvars[3+4*i+1]; /*2nd edge*/
    eField += ds[j] * ds[k] * cellvars[3+4*i+2]; /*3rd edge*/
    eField += (lcell[j]-ds[j]) * ds[k] * cellvars[3+4*i+3]; /*4th edge*/
    eField /= lcell[j]*lcell[k]; /*Divide field value over total face area*/

    bField = 0;
    bField += (lcell[i]-ds[i]) * cellvars[3+12+2*i]; /*1st face*/
    bField += ds[i] * cellvars[3+12+2*i+1]; /*2nd face*/
    bField /= lcell[i];

    if (i == 0) {
      field->e1 = eField;
      field->b1 = bField;
    } else if (i == 1) {
      field->e2 = eField;
      field->b2 = bField;
    } else {
      field->e3 = eField;
      field->b3 = bField;
    }
  }
}
