#ifndef INTERPOLATE_HPP
#define INTERPOLATE_HPP

#include "particle.hpp"

long n_cells; /* Number of grid cells */
double l_cell; /*Length of cell (edge length of cube)*/
double* cell_fields; /*Vector describing position of and all field elements of a cell*/
                     /*["least" corner vertex, E-field on edges, B-field on surfaces]*/
                     /*21 elements ordered as: [x1,x2,x3,E1,E2,E3,B1,B2,B3]*/
                     /*              of sizes:   1, 1, 1, 4, 4, 4, 2, 2, 2*/

int run_interpolation();
int get_gridvec_to_interpolate();
int interpolate_fields(double* pos, double* cfields, Field_part* fields); /*interpolate fields at particle position*/

#endif
