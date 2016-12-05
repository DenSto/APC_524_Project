#include <stdlib.h>
#include "grid.hpp"

Grid::Grid(int nxyz[3], int nGhosts, double xyz0[3], double dx): 
    nx_(nxyz[0]), 
    ny_(nxyz[1]), 
    nz_(nxyz[2]), 
    nGhosts_(nGhosts), 
    x0_(xyz0[0]), 
    y0_(xyz0[1]), 
    z0_(xyz0[2]), 
    dx_(dx), 
    iBeg_(nGhosts), 
    jBeg_(nGhosts), 
    kBeg_(nGhosts), 
    iEnd_(nx_-(nGhosts+1)), 
    jEnd_(ny_-(nGhosts+1)), 
    kEnd_(nz_-(nGhosts+1)),
    ghostVecSize_(9*(ny_-2*nGhosts)*(nz_-2*nGhosts))
{ 
    
    newField_(&Ex_); 
    newField_(&Ey_); 
    newField_(&Ez_); 
    newField_(&Bx_); 
    newField_(&By_); 
    newField_(&Bz_); 
    newField_(&Bx_tm1_); 
    newField_(&By_tm1_); 
    newField_(&Bz_tm1_); 
    newField_(&Jx_); 
    newField_(&Jy_); 
    newField_(&Jz_); 

    sliceTmp_ = new double[ghostVecSize_/9]; 
} 

Grid::~Grid() { 
    deleteField_(&Ex_); 
    deleteField_(&Ey_); 
    deleteField_(&Ez_); 
    deleteField_(&Bx_); 
    deleteField_(&By_); 
    deleteField_(&Bz_); 
    deleteField_(&Bx_tm1_); 
    deleteField_(&By_tm1_); 
    deleteField_(&Bz_tm1_); 
    deleteField_(&Jx_); 
    deleteField_(&Jy_); 
    deleteField_(&Jz_); 
    delete [] sliceTmp_; 
};

/* frees contiguous block fo memory for 3D array 
 * of size [nx_][ny_][nz_]. */ 
void Grid::newField_(double**** fieldPt) { 
    int i,j,k; // iterators 
    *fieldPt = new double**[nx_]; 
    for (i=0; i<nx_; ++i) { 
        *fieldPt[i] = new double*[ny_]; 
        for (j=0; j<ny_; ++j) { 
            *fieldPt[i][j] = new double[nz_]; 
        } 
    } 
};

/* frees contiguous block fo memory for 3D array 
 * of size [nx_][ny_][nz_]. */ 
void Grid::deleteField_(double**** fieldPt) { 
    int i,j,k; // iterators 
    for (i=0; i<nx_; ++i) { 
        for (j=0; j<ny_; ++j) { 
            delete [] *fieldPt[i][j]; 
        } 
        delete [] *fieldPt[i];
    } 
    delete [] *fieldPt;
}; 
