#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

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
#if USE_MPI
        double *getGlobalnyxz(void); // Global grid size
        double *getGlobalxyz0(void); // Global domain position
        double *getGlobalLxyz(void); // Global domain size
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
