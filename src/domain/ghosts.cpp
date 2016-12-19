#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "domain.hpp"

//! Allocate ghost buffers for MPI
/*! xgsize : size of ghost buffer in x direction
 *  ygsize : size of ghost buffer in y direction
 *  zgsize : size of ghost buffer in z direction
 */
void Domain::mallocGhosts(int xgsize, int ygsize, int zgsize){
    assert(xgsize>0);
    xghost_send_ = new double[2*xgsize]; // 1st half for left, 2nd for right
    xghost_recv_ = new double[2*xgsize]; // 1st half for left, 2nd for right
    assert(xghost_send_!=NULL);
    assert(xghost_recv_!=NULL);

    assert(ygsize>0);
    yghost_send_ = new double[2*ygsize]; // 1st half for left, 2nd for right
    yghost_recv_ = new double[2*ygsize]; // 1st half for left, 2nd for right
    assert(yghost_send_!=NULL);
    assert(yghost_recv_!=NULL);

    assert(zgsize>0);
    zghost_send_ = new double[2*zgsize]; // 1st half for left, 2nd for right
    zghost_recv_ = new double[2*zgsize]; // 1st half for left, 2nd for right
    assert(zghost_send_!=NULL);
    assert(zghost_recv_!=NULL);
};

void Domain::freeGhosts(void){
    delete[] xghost_send_;
    delete[] xghost_recv_;

    delete[] yghost_send_;
    delete[] yghost_recv_;

    delete[] zghost_send_;
    delete[] zghost_recv_;
};
