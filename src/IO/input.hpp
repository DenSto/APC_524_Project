#ifndef INPUT_HPP
#define INPUT_HPP

#if USE_MPI
#include "mpi.h"
#endif

//! Structure storing info in the input file
typedef struct {
    int nCell[3]; /*! number of cells in each direction */
    int nProc[3]; /*! number of processors to use in each direction */

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

    int nspecies; /// How many species of particles
                  /// eg. nspecies=2 in electron-proton plasma 

    long np; /// number of particles in each domain

    double t0;   /// start time of simulation

    double *mass_ratio;/// mass of each type of particle in unit of electron mass
                       /// array of length nspecies 
                       /// eg. in electron-proton plasma
                       ///     mass_ratio[0]=1; mass_ratio[1]=1830;

    double *charge_ratio;/// charge of each type of particle in unit of |e|
                         /// array of length nspecies
                         /// eq. in electron-proton plasma
                         ///   chargeratio[0]=-1; chargeratio[1]=1

    double *dens_frac; /// fractional density, array of length nspecies
                       /// eg. in quasineutral electron-proton plasma
                       ///     frac_dens[0]=0.5;
                       ///     frac_dens[1]=0.5;

    
    double *temp; /// Maxwellian temperature in unit of eV if specified
                  /// array of length nspecies
                  /// eq. in cold ion and hot electron plasma, possible value
                  ///     temp[0]=100; temp[1]=1.2;


    double xyz0[3];
    double Lxyz[3]; 

    // MPI can only send c_str!
    char distname[50]; /// name of file containing distribution function 
    char parts_bound[6][32]; /// particle boundary conditions for 6 sides of box
    char fields_bound[6][32];/// field boundary conditions for 6 sides of the box

} Input_Info_t;


/*! Class handeling input information */
class Input{
    public:
        Input(void);
        ~Input(void);

        int readinfo(char *fname);
        int checkinfo(void);
        void mallocinfo(int nspecies);
        Input_Info_t* getinfo(void);
#if USE_MPI
        void passinfo(int nspecies); 
#endif

    private:
        Input_Info_t *input_info_;
};

/*! Stand alone function to check input information passed by MPI */
void checkinput(Input_Info_t *input_info);


#endif
