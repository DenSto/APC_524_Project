#include <stdlib.h>
#include <math.h>
#include "interpolate.hpp"
#include "particle.hpp"

int get_gridvec_to_interpolate() {

}

int run_interpolation() {
  
}

int interpolate_fields(double* pos, double* cfields, Field_part* fields) {
  int i,j,k;
  double ds[3]; /*projected distances from particles to "least" cell faces*/
  double eField, bField; /*field variables for calculation*/

  for (i=0; i < 3; i++) {
    ds[i] = pos[i] - cfields[i];
  }

  for (i=0; i < 3; i++) {
    j = i+1 % 3; /*Modular arithmetic cycles over Cartesian directions*/
    k = i+2 % 3;

    eField = 0; /*Field additions are weighted by catercorner area.*/
    eField += (l_cell-ds[j]) * (l_cell-ds[k]) * cfields[3+4*i]; /*1st edge*/
    eField += ds[j] * (l_cell-ds[k]) * cfields[3+4*i+1]; /*2nd edge*/
    eField += ds[j] * ds[k] * cfields[3+4*i+2]; /*3rd edge*/
    eField += (l_cell-ds[j]) * ds[k] * cfields[3+4*i+3]; /*4th edge*/
    eField /= pow(l_cell,2.0); /*Divide field value over total face area*/

    bField = 0;
    bField += (l_cell-ds[i]) * cfields[3+12+2*i]; /*1st face*/
    bField += ds[i] * cfields[3+12+2*i+1]; /*2nd face*/
    bField /= l_cell;

    if (i == 0) {
      fields->e1 = eField;
      fields->b1 = bField;
    } else if (i == 1) {
      fields->e2 = eField;
      fields->b2 = bField;
    } else {
      fields->e3 = eField;
      fields->b3 = bField;
    }
  }
}
