#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "domain.hpp"

void Domain::mallocGhosts(Grid *grids){
    int xgsize = grids->getGhostVecSize(); // ghost size in x direction
    assert(xgsize>0);
    xghost_send_ = (double*)malloc(2*xgsize); // 1st half for left, 2nd for right
    xghost_recv_ = (double*)malloc(2*xgsize); // 1st half for left, 2nd for right
    assert(xghost_send_!=NULL);
    assert(xghost_recv_!=NULL);
};

void Domain::freeGhosts(void){
    free(xghost_send_);
    free(xghost_recv_);
};
