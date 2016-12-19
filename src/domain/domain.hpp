#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

#include "../grid/grid.hpp"
#include "../IO/IO.hpp"

#if USE_MPI
#include "mpi.h"
#endif

class Domain {
    public:
        Domain(Input_Info_t *input_info);
        ~Domain();
        
        int getnGhosts(void);
        int *getnxyz(void);    // Local grid size
        int *getn2xyz(void);    // Local grid area
        double *getxyz0(void); // Local origin
        double *getLxyz(void); // Local domain size
        double getmindx(void);

        void mallocGhosts(int xgsize, int ygsize, int zgsize); // allocate ghostVec
        void freeGhosts(void); // free ghostVec

        void PassFields(Grid *grids, Input_Info_t *input_info); // field boundary
 
	int *getnProcxyz(void); // return pointer nProcxyz_
        int *getmyijk(void); 
		int *getNeighbours(); // like below, but in useful array form xLR-yLR-zLR
        int getxl(void); // rank of x left neighbot 
        int getyl(void); 
        int getzl(void);
        int getxr(void); // rank of x right neighbor
        int getyr(void);
        int getzr(void);

        // return rank for assigned i,j,k
	int ijkToRank(int i, int j, int k); 
        // assign value to allocated myijk[3] 
        void RankToijk(int rank, int *myijk); 
	
    private:
        int size_; // MPI size
        int rank_; // MPI rank     

        int nGhosts_; 
        int *nxyz_; // nxyz_[0]=nx
        int *n2xyz_; // n2xyz_[0]=ny*nz
        int n3xyz_; // n3xyz=nx*ny*nz
        double *xyz0_;
        double *Lxyz_;

        int rank_xl_, rank_xr_; // x neighbors 
        int rank_yl_, rank_yr_; // y neighbors
        int rank_zl_, rank_zr_; // z neighbors

        double *xghost_send_, *xghost_recv_; // buffer for ghostVec in x direction
        double *yghost_send_, *yghost_recv_; // buffer for ghostVec in y direction
        double *zghost_send_, *zghost_recv_; // buffer for ghostVec in z direction

        int *nProcxyz_;   // sizes of the partitions
		int *neighbours_;
		int *myijk_; // where am I on this partition?

#if USE_MPI
        // communication in x direction
        MPI_Request xreqs_[4];
        MPI_Status xstats_[4];
#endif
};

void checkdomain(Domain *domain);

#endif
