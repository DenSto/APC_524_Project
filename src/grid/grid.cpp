#include <stdio.h> 
#include <stdlib.h>
#include <assert.h> 
#include "grid.hpp"
#include "../IO/IO.hpp"

/// Grid constructor 
/*! Input arguments: \n 
 * nxyz: integer array [nx,ny,nz] where nx is the total number of cells (physical + ghost) in the x direction in the simulation, and the same for ny,nz. \n 
 * nGhosts: integer number of ghost cells on each side of the domain. This should always be at least 1. Currently the code does not support nGhosts>1, though it may in the future (to take advantage of higher order finite difference and interpolation methods, for instance). \n 
 * xyz0: integer array [x0,y0,z0] where x0 is the initial x position, and the same for y0,z0 \n 
 * Lxyz0: double array [Lx,Ly,Lz] where Lx is the physical length of each cell in the x direction, and the same for Ly,Lz \n 
 */ 
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
    dx_(Lxyz[0]/nxyz[0]), 
    dy_(Lxyz[1]/nxyz[1]), 
    dz_(Lxyz[2]/nxyz[2]),
    idx_(1.0/dx_), 
    idy_(1.0/dy_),
    idz_(1.0/dz_),
    nRealPtsYZPlane_((ny_+1-2*nGhosts)*(nz_+1-2*nGhosts)), // fields have ni_+1-2*nGhosts physical points in ith direction
    nFields_(9),
    ghostVecSize_(nFields_*nRealPtsYZPlane_)
{
    checkInput_(); 
    
    Ex_=newField_(); 
    Ey_=newField_(); 
    Ez_=newField_(); 
    Bx_=newField_(); 
    By_=newField_(); 
    Bz_=newField_(); 
    Bx_tm1_=newField_(); 
    By_tm1_=newField_(); 
    Bz_tm1_=newField_(); 
    rhox_=newField_(); 
    rhoy_=newField_(); 
    rhoz_=newField_(); 
    Jx_=newField_(); 
    Jy_=newField_(); 
    Jz_=newField_(); 

    sliceTmp_ = new double[ghostVecSize_/nFields_]; 
} 

/// Grid destructor 
/*! calls deleteField_ on each of the double*** fields 
 */ 
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
    deleteField_(rhox_); 
    deleteField_(rhoy_); 
    deleteField_(rhoz_); 
    deleteField_(Jx_); 
    deleteField_(Jy_); 
    deleteField_(Jz_); 

    delete [] sliceTmp_; 
};

/// allocates contiguous block of memory for a single field 
/*! Returns double*** of size [nx_+1][ny_+1][nz_+1]
 */ 
double*** Grid::newField_() { 
    int i,j; // iterators 
    double*** fieldPt = new double**[nx_+1]; 
    assert( fieldPt != NULL ); 
    for (i=0; i<nx_+1; ++i) { 
        fieldPt[i] = new double*[ny_+1]; 
        assert( fieldPt[i] != NULL ); 
        for (j=0; j<ny_+1; ++j) { 
            fieldPt[i][j] = new double[nz_+1]; 
            assert( fieldPt[i][j] != NULL ); 
        } 
    } 
    return fieldPt;
};

/// frees contiguous block of memory for a single field
/*! Deletes double*** of size [nx_+1][ny_+1][nz_+1]
*/ 
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

/// checks validity of input parameters for Grid constructor 
/*! asserts necessary conditions on each input (mainly positivity of many parameters). Terminates program if inputs are incorrect.
 */ 
void Grid::checkInput_() { 
    assert(nx_ > 0); 
    assert(ny_ > 0); 
    assert(nz_ > 0); 
    assert(nGhosts_ == 1); // currently some grid functions assume this, though they can be generalized later to allow fo rmore 
    assert(Lx_ > 0); 
    assert(Ly_ > 0); 
    assert(Lz_ > 0); 
    assert(iBeg_ > 0); 
    assert(jBeg_ > 0); 
    assert(kBeg_ > 0); 
    assert(iEnd_ < nx_+1); 
    assert(jEnd_ < ny_+1); 
    assert(kEnd_ < nz_+1);
    assert(dx_ > 0); 
    assert(dy_ > 0); 
    assert(dz_ > 0);
    assert(nRealPtsYZPlane_ > 0); 
    assert(nFields_ == 9); 
    assert(ghostVecSize_ > 0); 
}; 

/// sets all of J (Jx,Jy,Jz) to be identically zero
void Grid::zeroJ() { 
    int i,j,k; // iterators 
    for (i=0; i<nx_+1; ++i) { 
        for (j=0; j<ny_+1; ++j) { 
            for (k=0; k<nz_+1; ++k) { 
                Jx_[i][j][k]=0; 
                Jy_[i][j][k]=0; 
                Jz_[i][j][k]=0;
            } 
        } 
    } 
} 

/// Initialize E and B fields
/*! Use restart file to set values of initial E,B,J fields
 */ 
void Grid::InitializeFields(int restart){
}; 
