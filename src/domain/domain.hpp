#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

#include "../grid/grid.hpp"
#include "../particles/particle_handler.hpp"
#include "../IO/IO.hpp"

class Domain {
    public:
        Domain(int size, int rank, Input_Info_t *input_info);
        ~Domain();
        
        int getnGhosts(void);
        int *getnxyz(void);    // Local grid size
        double *getxyz0(void); // Local position
        double *getLxyz(void); // Local domain size
        double getmindx(void);

        void mallocGhosts(void); // allocate ghostVec
        void freeGhosts(void); // free ghostVec

        void PassParticles(Particle_Handler *parts_fields); // particle boundary
        void PassFields(Grid *grids); // field boundary
 
#if USE_MPI
        int *getGlobalnyxz(void); // Global grid size
        double *getGlobalxyz0(void); // Global domain position
        double *getGlobalLxyz(void); // Global domain size
		int *getMyLocationOnMap();
		int *getnProcxyz(void);
		int ***getProcMap(void);
#endif
	
    private:
        int size_; // MPI size
        int rank_; // MPI rank     

        int nGhosts_; 
        int *nxyz_;
        double *xyz0_;
        double *Lxyz_;

        int nxny_; //size of z-plane nx*ny
        int nynz_; //size of x-plane ny*nz
        int nznx_; //size of y-plane nz*nx

        int rank_xl_, rank_xr_; // x neighbors 
        int rank_yl_, rank_yr_; // y neighbors
        int rank_zl_, rank_zr_; // z neighbors

        double *xghost_send, *xghost_recv; // buffer for ghostVec in x direction
        double *yghost_send, *yghost_recv; // buffer for ghostVec in y direction
        double *zghost_send, *zghost_recv; // buffer for ghostVec in z direction


#if USE_MPI
		int ***procMap_;  // map of rank to partition segment
		int *nProcxyz_;   // sizes of the partitions
		int *myLocationOnMap_; // where am I on this partition?
        int *globalnxyz_;
        double *globalxyz0_;
        double *globalLxyz_;
#endif
};

void checkdomain(int rank, Domain *domain);

#endif
