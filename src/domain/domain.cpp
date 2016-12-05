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
       Lxyz_ new double[3];
       Lxyz_[0] = 1.0/size;
       Lxyz_[1] = 1.0;
       Lxyz_[2] = 1.0;
  
}

Domain::~Domain(){
    delete[] nxyz_;
    delete[] xyz0_;
    delete[] Lxyz_;
}
