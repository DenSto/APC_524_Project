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

    private:
        int size_; // MPI size
        int rank_; // MPI rank
        int nGhosts_; 
        int *nxyz_;
        double *xyz0_;
        double *Lxyz_;
};

void checkdomain(int rank, Domain domain);

#endif
