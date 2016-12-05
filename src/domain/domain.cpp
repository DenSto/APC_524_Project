#include "domain.hpp"

Domain::Domain(int size, int rank, Input_Info_t *input_info)
      : size_(size),
        rank_(rank){

       nGhosts_ = 1;

       const int nx = input_info->nx; 
       nxyz_ = new int[3];
       nxyz_[0]=nx/size;
       nxyz_[1]=nx;
       nxyz_[2]=nx;

       //const double *xyz = input_info->xyz; 
       xyz0_ = new double[3];
       xyz0_[0] = 1.0*rank;
       xyz0_[1] = 0.0;
       xyz0_[2] = 0.0;

       //const double *lxyz = input_info->lxyz;
       Lxyz_ = new double[3];
       Lxyz_[0] = 1.0/size;
       Lxyz_[1] = 1.0;
       Lxyz_[2] = 1.0;
  
}

Domain::~Domain(){
    delete[] nxyz_;
    delete[] xyz0_;
    delete[] Lxyz_;
}

int Domain::getnGhosts(void){
    return nGhosts_;
}

int *Domain::getnxyz(void){
    return nxyz_;
}

double *Domain::getxyz0(void){
    return xyz0_;
}

double *Domain::getLxyz(void){
    return Lxyz_;
}

void checkdomain(int rank, Domain domain){
    printf("rank=%d,nGhosts=%d\n",rank,domain.getnGhosts());
    int *nxyz = domain.getnxyz();
    printf("rank=%d,nxyz=%d,%d,%d\n",rank,nxyz[0],nxyz[1],nxyz[2]);

    double *xyz0 = domain.getxyz0();
    printf("rank=%d,xyz0=%f,%f,%f\n",rank,xyz0[0],xyz0[1],xyz0[2]);

    double *Lxyz = domain.getLxyz();
    printf("rank=%d,Lxyz=%f,%f,%f\n",rank,Lxyz[0],Lxyz[1],Lxyz[2]);
}
    
