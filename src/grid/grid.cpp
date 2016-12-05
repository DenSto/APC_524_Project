#include <stdlib.h>
#include "grid.hpp"

Grid::Grid(int *nxyz, int nGhosts, double *xyz0, double *Lxyz): 
    nx_(nxyz[0]), 
    ny_(nxyz[1]), 
    nz_(nxyz[2]), 
    nGhosts_(nGhosts), 
    x0_(xyz0[0]), 
    y0_(xyz0[1]), 
    z0_(xyz0[2]), 
    Lx_(Lxyz[0]), 
    Ly_(Lxyz[1]), 
    Lz_(Lxyz[2]), 
    iBeg_(nGhosts), 
    jBeg_(nGhosts), 
    kBeg_(nGhosts), 
    iEnd_(nx_-nGhosts), // fields are length nx_+1, so last element is indexed as nx_, then subtract nGhosts to get index of last physical point
    jEnd_(ny_-nGhosts), 
    kEnd_(nz_-nGhosts),
    nFields_(9),
    nRealPtsYZPlane_((ny_+1-2*nGhosts)*(nz_+1-2*nGhosts)),
    dx_(Lxyz[0]/nxyz[0]), 
    dy_(Lxyz[1]/nxyz[1]), 
    dz_(Lxyz[2]/nxyz[2]),
    // fields have ni_+1-2*nGhosts physical points in ith direction
    ghostVecSize_(nFields_*nRealPtsYZPlane_)
{ 
    
    Ex_=newField_(); 
    Ey_=newField_(); 
    Ez_=newField_(); 
    Bx_=newField_(); 
    By_=newField_(); 
    Bz_=newField_(); 
    Bx_=newField_(); 
    By_=newField_(); 
    Bz_=newField_(); 
    Jx_=newField_(); 
    Jy_=newField_(); 
    Jz_=newField_(); 

    sliceTmp_ = new double[ghostVecSize_/nFields_]; 
} 

Grid::~Grid() { 
    deleteField_(Ex_); 
    deleteField_(Ey_); 
    deleteField_(Ez_); 
    deleteField_(Bx_); 
    deleteField_(By_); 
    deleteField_(Bz_); 
    deleteField_(Bx_tm1_); 
    deleteField_(By_tm1_); 
    deleteField_(Bz_tm1_); 
    deleteField_(Jx_); 
    deleteField_(Jy_); 
    deleteField_(Jz_); 

    delete [] sliceTmp_; 
};

/* allocates contiguous block of memory for 3D array 
 * of size [nx_][ny_][nz_]. */ 
double*** Grid::newField_() { 
    int i,j; // iterators 
    double*** fieldPt = new double**[nx_+1]; 
    for (i=0; i<nx_+1; ++i) { 
        fieldPt[i] = new double*[ny_+1]; 
        for (j=0; j<ny_+1; ++j) { 
            fieldPt[i][j] = new double[nz_+1]; 
        } 
    } 
    return fieldPt;
};

/* frees contiguous block fo memory for 3D array 
 * of size [nx_][ny_][nz_]. */ 
void Grid::deleteField_(double*** fieldPt) { 
    int i,j; // iterators 
    for (i=0; i<nx_+1; ++i) { 
        for (j=0; j<ny_+1; ++j) { 
            delete [] fieldPt[i][j]; 
        } 
        delete [] fieldPt[i];
    } 
    delete [] fieldPt;
}; 
