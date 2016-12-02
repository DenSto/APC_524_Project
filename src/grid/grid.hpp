#ifndef GRID_HPP
#define GRID_HPP

typedef struct {
    double *Efield;
    double *Bfield;
} Field_grid_t;

typedef struct {
    double *Rhofield;
    double *Jfield;
} RhoJ_grid_t;

typedef struct {
    Field_grid_t Fields;
    RhoJ_grid_t  RhoJ;
} Grids_list_t;

#endif
