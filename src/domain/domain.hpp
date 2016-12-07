#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

#include "../IO/IO.hpp"

class Domain {
    public:
        Domain(int size, int rank, Input_Info_t *input_info);
        ~Domain();
        
        int getnGhosts(void);
        int *getnxyz(void);
        double *getxyz0(void);
        double *getLxyz(void);
        double getmindx(void); 
#if USE_MPI
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
		int *myLocationOnMap; // where am I on this partition?
#endif
};

void checkdomain(int rank, Domain *domain);

#endif
