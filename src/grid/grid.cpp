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
    nxTot_(nx_+1), 
    nyTot_(ny_+1), 
    nzTot_(nz_+1),
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
    nRealPtsYZPlane_((nyTot_-2*nGhosts)*(nzTot_-2*nGhosts)), // fields have ni_+1-2*nGhosts physical points in ith direction
    nFieldsToSend_(9),
    nFieldsTotal_(15), 
    ghostVecSize_(nFieldsToSend_*nRealPtsYZPlane_)
{
    checkInput_(); 

    sliceTmp_ = new double[ghostVecSize_/nFieldsToSend_]; 
    fieldIsContiguous_ = new double[nFieldsTotal_];
 
    int ifield = -1; 
    Ex_=newField_(++ifield); 
    Ey_=newField_(++ifield); 
    Ez_=newField_(++ifield); 
    Bx_=newField_(++ifield); 
    By_=newField_(++ifield); 
    Bz_=newField_(++ifield); 
    Jx_=newField_(++ifield); 
    Jy_=newField_(++ifield); 
    Jz_=newField_(++ifield); 
    Bx_tm1_=newField_(++ifield); 
    By_tm1_=newField_(++ifield); 
    Bz_tm1_=newField_(++ifield); 
    rhox_=newField_(++ifield); 
    rhoy_=newField_(++ifield); 
    rhoz_=newField_(++ifield); 
} 

/// Grid destructor 
/*! calls deleteField_ on each of the double*** fields 
 */ 
Grid::~Grid() { 
    /* note: these must be deleted in the same order as they were created
     * since they use fieldIsContiguous_ to determine the create deletion method (contiguous vs noncontiguous) */ 
    int ifield=-1; 
    deleteField_(Ex_,++ifield); 
    deleteField_(Ey_,++ifield); 
    deleteField_(Ez_,++ifield); 
    deleteField_(Bx_,++ifield); 
    deleteField_(By_,++ifield); 
    deleteField_(Bz_,++ifield); 
    deleteField_(Jx_,++ifield); 
    deleteField_(Jy_,++ifield); 
    deleteField_(Jz_,++ifield); 
    deleteField_(Bx_tm1_,++ifield); 
    deleteField_(By_tm1_,++ifield); 
    deleteField_(Bz_tm1_,++ifield); 
    deleteField_(rhox_,++ifield); 
    deleteField_(rhoy_,++ifield); 
    deleteField_(rhoz_,++ifield); 

    delete [] sliceTmp_;
    delete [] fieldIsContiguous_; 
};

/// allocates memory for a single field 
/*! Returns double*** of size [nx_+1][ny_+1][nz_+1]. \n
 * First attempts to allocate contiguously. If that fails, issues a warning and attempts to allocate with several calls to new. 
 */ 
double*** Grid::newField_(int ifield) { 
    int i,j; // iterators 
    
    // try to allocate a contiguous block in memory 
    double*** fieldPt = new double** [nxTot_]; 
    fieldPt[0] = new double* [nxTot_*nyTot_]; 
    fieldPt[0][0] = new double [nxTot_*nyTot_*nzTot_]; 
    if (fieldPt != NULL && fieldPt[0] != NULL && fieldPt[0][0] != NULL) { 
        fieldIsContiguous_[ifield] = 1; 
        for (i=0; i<nxTot_; ++i) { 
            if (i < nxTot_-1) { 
                fieldPt[0][(i+1)*nyTot_] = &(fieldPt[0][0][(i+1)*nyTot_*nzTot_]); 
                fieldPt[i+1] = &(fieldPt[0][(i+1)*nyTot_]); 
            } 
            for (j=0; j<nyTot_; ++j) { 
                if (j > 0) { 
                    fieldPt[i][j] = fieldPt[i][j-1] + nzTot_; 
                } 
            } 
        } 
    } 
    else { 
        // if contiguous allocation failed, use a non contiguous allocation
        // e.g. there could be enough fragmented memory to use
        printf("WARNING: unable to allocate fields contiguously. Attempting to allocate noncontiguously. This may impact performance."); 
        fieldIsContiguous_[ifield] = 0; 
        double*** fieldPt = new double**[nxTot_]; 
        assert( fieldPt != NULL ); 
        for (i=0; i<nxTot_; ++i) { 
            fieldPt[i] = new double*[nyTot_]; 
            assert( fieldPt[i] != NULL ); 
            for (j=0; j<nyTot_; ++j) { 
                fieldPt[i][j] = new double[nzTot_]; 
                assert( fieldPt[i][j] != NULL ); 
            } 
        } 
    } 
    return fieldPt; 
};

/// frees memory for a single field
/*! Uses fieldIsContiguous_ to determine contiguous or noncontiguous deltion method
*/ 
void Grid::deleteField_(double*** fieldPt, int ifield) { 
    int i,j; // iterators 
    if (fieldIsContiguous_[ifield] == 1) { 
        delete [] fieldPt[0][0]; 
        delete [] fieldPt[0]; 
        delete [] fieldPt; 
    } 
    else { 
        // use deletion method corresponding to the noncontiguous memory allocation method 
        for (i=0; i<nxTot_; ++i) { 
            for (j=0; j<nyTot_; ++j) { 
                delete [] fieldPt[i][j]; 
            } 
            delete [] fieldPt[i];
        } 
        delete [] fieldPt;
    } 
};

/// checks validity of input parameters for Grid constructor 
/*! asserts necessary conditions on each input (mainly positivity of many parameters). Terminates program if inputs are incorrect.
 */ 
void Grid::checkInput_() { 
    assert(nx_ > 0); 
    assert(ny_ > 0); 
    assert(nz_ > 0); 
    assert(nGhosts_ == 1); // currently some grid functions assume this, though they can be generalized later to allow fo rmore 
    assert(nxTot_ > nx_); 
    assert(nyTot_ > ny_); 
    assert(nzTot_ > nz_); 
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
    assert(nFieldsToSend_ == 9); 
    assert(nFieldsTotal_ == 15); 
    assert(ghostVecSize_ > 0); 
}; 

/// sets all of J (Jx,Jy,Jz) to be identically zero
/*! Used during particle deposition. 
 */ 
void Grid::zeroJ() { 
    int i,j,k; // iterators 
    for (i=0; i<nxTot_; ++i) { 
        for (j=0; j<nyTot_; ++j) { 
            for (k=0; k<nzTot_; ++k) { 
                Jx_[i][j][k]=0; 
                Jy_[i][j][k]=0; 
                Jz_[i][j][k]=0;
            } 
        } 
    } 
}; 

/// Initialize E and B fields
/*! Use restart file to set values of initial E,B,J fields
 */ 
void Grid::InitializeFields(int restart){
    /* jlestz: dummy code setting all fields to zero
     * placeholder until restart file format is determined */ 
    int i,j,k; 
    for (i=0; i<nxTot_; ++i) { 
        for (j=0; j<nyTot_; ++j) { 
            for (k=0; k<nzTot_; ++k) { 
                Ex_[i][j][k]=0; 
                Ey_[i][j][k]=0; 
                Ez_[i][j][k]=0; 
                Bx_[i][j][k]=0; 
                By_[i][j][k]=0; 
                Bz_[i][j][k]=0; 
                Jx_[i][j][k]=0; 
                Jy_[i][j][k]=0; 
                Jz_[i][j][k]=0; 
            } 
        } 
    } 
}; 
