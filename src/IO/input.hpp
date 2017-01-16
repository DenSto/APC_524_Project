#ifndef INPUT_HPP
#define INPUT_HPP

#if USE_MPI
#include "mpi.h"
#endif

#define NDIM 3 // spatial dimension
#define NCHAR 64 // characters reserved for MPI
#define NSPEC 32 // max nspecies reserved for MPI_BCast
                 // MPI cannot send structure containing unallocated pointers
#define NWAVE 32 // max nwaves for MPI_BCase

//! Structure storing info in the input file
typedef struct {
    int nCell[NDIM]; /*! number of cells in each direction */
    int nProc[NDIM]; /*! number of processors to use in each direction */

    int nt; /*! number of time steps */
    int restart; /// How many previous runs?
                 /// restart = 0: initial run
                 //  restart = 3: third restart run
    int debug; /// 0: do not print debug statements
               /// 1: print minimal debug statements
               /// 2: print more debug statements
               /// 3: write debug files
    int relativity; /// 1: use relativistic pusher
                    /// 0: use nonrelativistic pusher
    int electrostatic; /// 1: use electrostatic field solve
                       /// 0: use electromagnetic field solve 

    int nspecies; /// How many species of particles
                  /// eg. nspecies=2 in electron-proton plasma 
                  /// nspecies <=NSPEC

    // diagnostics parameters
    // number of time steps between writing files    
    int nstep_fields;  
    int nstep_parts;  
    int nstep_restart;
    int nstep_sort;

    int which_fields; /// flag determine which fields to write
                      /// 0: write components of rho
                      /// 1: write components of E 
                      /// 2: write components of B
                      /// 3: write components of J
                      /// 4: write all fields
    int output_pCount; /// how many particles per core to print

    int nwaves; /// How many waves to inject into the system
                /// nwave<=NWAVE
    int inSide[NWAVE]; /// from which sides are waves injected
                       /// eg. inSide[0]=-1: 1st wave injected in x direction(1) from left(-)
                       ///     inside[1]=+3: 2nd wave injected in z direction(3) from right(+)
    int inPolE[NWAVE]; /// polarization of E field of injected waves
                       /// eg. inPolE[0]=2: 1st wave E field is in y direction(2)
                       ///     inPolE[1]=3: 2nd wave E field is in Z direction(3)
                       ///     inPolE should only take value of 1,2,3

    int isTestParticle[NSPEC];/// is the species a test particle species
                       /// i.e. it feels fields but does not influence them
                       /// 0 for no, 1 for yes
    long nparticles_tot; /// total number of particles of all species in the entire 
                         /// simulation box
    long nparticles_domain; /// total number of particles of all species in each domain

    double t0;   /// start time of simulation

    double dens_phys; /// physical number density of all particles
                      /// used to scale mass, charge and temperature
                      /// of super particles.

    double super_ratio; /// ratio of physical density over PIC density
                        /// super_ratio = -1 when there is not particle in simulation

    double mass_ratio[NSPEC];/// mass of each type of particle in unit of electron mass
                       /// array of length nspecies 
                       /// eg. in electron-proton plasma
                       ///     mass_ratio[0]=1; mass_ratio[1]=1830;
                       //

    double charge_ratio[NSPEC];/// charge of each type of particle in unit of |e|
                         /// array of length nspecies
                         /// eq. in electron-proton plasma
                         ///   chargeratio[0]=-1; chargeratio[1]=1

    double dens_frac[NSPEC]; /// fractional density, array of length nspecies
                       /// eg. in quasineutral electron-proton plasma
                       ///     frac_dens[0]=0.5;
                       ///     frac_dens[1]=0.5;

    
    double temp[NSPEC]; /// Maxwellian temperature in unit of eV if specified
                  /// array of length nspecies
                  /// eq. in cold ion and hot electron plasma, possible value
                  ///     temp[0]=100; temp[1]=1.2;

    double peakamps[NWAVE]; // wave amplitudes(in program unit) injected at the boundary
    double omegas[NWAVE]; // corresponding wave frequencies (in program unit)
    double phases[NWAVE]; // corresponding wave phase (in unit of rad)

    double invWidths[NWAVE]; // the inverse of Gaussian pulse width in time (in program unit)
    double delays[NWAVE]; // time delay of Gaussian pulse center w.r.t. time_phys

    double E0[NDIM]; /// background electric field
    double B0[NDIM]; /// background magnetic field

    double bound_phi[2*NDIM]; /// boundary conditions for poisson initialization
    double bound_Ax[2*NDIM]; /// boundary conditions for poisson initialization
    double bound_Ay[2*NDIM]; /// boundary conditions for poisson initialization
    double bound_Az[2*NDIM]; /// boundary conditions for poisson initialization

    double xyz0[NDIM]; // origin of simulation
    double Lxyz[NDIM]; // physical length of simulation

    // MPI can only send c_str!
    char distname[NCHAR]; /// name of file containing distribution function 
    
    char parts_init[NCHAR]; /// particle initialization method
    char fields_init[NCHAR]; /// field initialization method

    char parts_bound[2*NDIM][NCHAR]; /// particle boundary conditions for 6 sides of box
    char fields_bound[2*NDIM][NCHAR];/// field boundary conditions for 6 sides of the box

} Input_Info_t;


/*! Class handeling input information */
class Input{
    public:
        Input(void);
        ~Input(void);

        int readinfo(char* inputname);
        int ProcessInfo(void);
        Input_Info_t* getinfo(void);
#if USE_MPI
        void passinfo(void); 
#endif

    private:
        Input_Info_t *input_info_;
};

/*! Stand alone function to check input information passed by MPI */
void checkinput(Input_Info_t *input_info);


#endif
