#ifndef RESOLVE_HPP_
#define RESOLVE_HPP_

#include "../globals.hpp"
#include "../IO/input.hpp"
#include "../grid/grid.hpp"
#include "../particles/particle_handler.hpp"

#if USE_MPI
#include "mpi.h"
#endif

//! Sturcture storing time resolutions information, in unit of ps or THz
typedef struct { 
    double time_light; //light transit unit cell
    double omega_p; //plasma frequency
    double omega_c; //maximum gyro frequency
    double omega_e; //maximum boundary wave frequency   
 
} Time_Resolve_t;

//! Sturcture storing spatial resolutions information, in unit of cm or cm^{-1}
typedef struct {
    double dCell; // minimum length of unit cells
    double Debye; // Debye length
    double Skindepth; // Skin depth;
    double wavelength; // minimum external vacuum EM wave length
} Space_Resolve_t;

//! Structure storing both time and space resolution
typedef struct {
    Time_Resolve_t time;
    Space_Resolve_t space;
} Resolve_t;

class Resolution {
    public:
        Resolution(Input_Info_t *input_info);
        ~Resolution(void);

        double DetermineTimeStep(void);
        double DetermineCellSize(void);
    private:
        Resolve_t *resolve_;
};

#endif
