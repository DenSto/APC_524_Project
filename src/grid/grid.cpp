#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <assert.h> 
#include "grid.hpp"
#include "../globals.hpp"
#include "../IO/input.hpp"

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
    dx_(Lxyz[0]/nxyz[0]), 
    dy_(Lxyz[1]/nxyz[1]), 
    dz_(Lxyz[2]/nxyz[2]),
    idx_(1.0/dx_), 
    idy_(1.0/dy_),
    idz_(1.0/dz_),
    maxPointsInPlane_(std::max(std::max(nxTot_*nyTot_,nxTot_*nzTot_),nyTot_*nzTot_)),
    nFieldsToSend_(9),
    ghostVecSize_(nFieldsToSend_*maxPointsInPlane_), 
    nFieldsTotal_(21), 
    ExID_(0),
    EyID_(1),
    EzID_(2),
    BxID_(3),
    ByID_(4),
    BzID_(5),
    JxID_(6),
    JyID_(7),
    JzID_(8),
    Bx_tm1ID_(9),
    By_tm1ID_(10),
    Bz_tm1ID_(11),
    rhoID_(12),
    nTypes_(7), 
    edgeXID_(0), 
    edgeYID_(1), 
    edgeZID_(2), 
    faceXID_(3),
    faceYID_(4),
    faceZID_(5),
    vertID_(6)
{
   
    checkInput_(); 

    fieldIsContiguous_ = new double[nFieldsTotal_];
 
    Ex_=newField_(ExID_); 
    Ey_=newField_(EyID_); 
    Ez_=newField_(EzID_); 
    Bx_=newField_(BxID_); 
    By_=newField_(ByID_); 
    Bz_=newField_(BzID_); 
    Jx_=newField_(JxID_); 
    Jy_=newField_(JyID_); 
    Jz_=newField_(JzID_); 
    Bx_tm1_=newField_(Bx_tm1ID_); 
    By_tm1_=newField_(By_tm1ID_); 
    Bz_tm1_=newField_(Bz_tm1ID_); 
    rho_=newField_(rhoID_); 
    
    fieldType_ = setFieldType_(); 
    fieldSize_ = setFieldSize_(); 
    fieldPtr_ = setFieldPtr_(); 

    sliceTmp = new double[maxPointsInPlane_]; 
    ghostTmp = new double[ghostVecSize_]; 
} 

/// Grid destructor 
/*! calls deleteField_ on each of the double*** fields 
 */ 
Grid::~Grid() { 
    deleteField_(Ex_,ExID_); 
    deleteField_(Ey_,EyID_); 
    deleteField_(Ez_,EzID_); 
    deleteField_(Bx_,BxID_); 
    deleteField_(By_,ByID_); 
    deleteField_(Bz_,BzID_); 
    deleteField_(Jx_,JxID_); 
    deleteField_(Jy_,JyID_); 
    deleteField_(Jz_,JzID_); 
    deleteField_(Bx_tm1_,Bx_tm1ID_); 
    deleteField_(By_tm1_,By_tm1ID_); 
    deleteField_(Bz_tm1_,Bz_tm1ID_); 
    deleteField_(rho_,rhoID_); 

    delete [] fieldIsContiguous_; 
    deleteFieldType_(); 
    deleteFieldSize_(); 
    deleteFieldPtr_();
    delete [] sliceTmp; 
    delete [] ghostTmp; 
};

/// allocates memory for a single field 
/*! Returns double*** of size [nx_+1][ny_+1][nz_+1]. \n
 * First attempts to allocate contiguously. If that fails, issues a warning and attempts to allocate with several calls to new. 
*/ 
double*** Grid::newField_(int fieldID) { 
    int i,j; // iterators 
    
    // try to allocate a contiguous block in memory 
    double*** fieldPt = new double** [nxTot_]; 
    fieldPt[0] = new double* [nxTot_*nyTot_]; 
    fieldPt[0][0] = new double [nxTot_*nyTot_*nzTot_]; 
    if (fieldPt != NULL && fieldPt[0] != NULL && fieldPt[0][0] != NULL) { 
        fieldIsContiguous_[fieldID] = 1; 
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
        fieldIsContiguous_[fieldID] = 0; 
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
void Grid::deleteField_(double*** fieldPt, int fieldID) { 
    int i,j; // iterators 
    if (fieldIsContiguous_[fieldID] == 1) { 
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

/// constructs and returns fieldSize_ array 
/*! fieldSize_ is an ntypes by ndim array storing the number 
 * of physical + ghost points in each direction. This is necessary 
 * because although all field arrays are allocated to be the same 
 * size (nx+1,ny+1,nz+1), due to the different locations of each 
 * type of field on the grid (3 types of edge locations, 3 types of 
 * face locations, vertices) which leads to differences in the number
 * of points needed for nx,ny,nz cells. \n
 * rows correspond to fieldType: 
 * 0: x edge (Ex/Jx), 1: y edge (Ey/Jy), 2: z edge (Ez/Jz), \n
 * 3: x face (Bx)   , 4: y face (By)   , 5: z face (Bz), \n
 * 6: vertices (rho) \n
 * columns correspond to the direction (0,1,2)=(x,y,z)
 */ 
int** Grid::setFieldSize_() { 
    int i,j; 
    int ndim=3; 
    int** fieldSize = new int*[nTypes_]; 
    assert(fieldSize != NULL); 
    for (i=0; i<nTypes_; ++i) { 
        fieldSize[i] = new int[ndim]; 
        assert(fieldSize[i] != NULL); 
    };

    // set the actual values 
    int nxyz[3] = {nx_,ny_,nz_}; 
    int edge,dir; 
    for (i=0; i<nTypes_; ++i) { 
        if (i < ndim) { 
            edge=1; 
        }
        else { 
            edge=0; 
        };
        dir = (i % ndim); 
        for (j=0; j<ndim; ++j) { 
            fieldSize[i][j] = nxyz[j]+edge; 
            if (i < nTypes_-1) { 
                if (j == dir) { 
                    fieldSize[i][j] += pow(-1,edge); 
                }; 
            } 
            else { 
                ++fieldSize[i][j]; 
            }; 
        }; 
    };
    return fieldSize; 
};

///  deletes fieldSize_ array 
void Grid::deleteFieldSize_() { 
    int i; 
    for (i=0; i<nTypes_; ++i) { 
        delete [] fieldSize_[i]; 
    }; 
    delete [] fieldSize_; 
}; 

/// constructs and returns fieldType_ array 
/*! fieldType_ is an nFieldsTotal_ array of ints storing 
 * the type of each field (edgeX, faceZ, vertex, etc). \n
 * e.g. int typeOfBx = fieldType_[BxID_]; 
 */ 
int* Grid::setFieldType_() { 
    int* fieldType = new int[nFieldsTotal_];
    
    fieldType[ExID_] = edgeXID_; 
    fieldType[EyID_] = edgeYID_; 
    fieldType[EzID_] = edgeZID_; 
    fieldType[BxID_] = faceXID_; 
    fieldType[ByID_] = faceYID_; 
    fieldType[BzID_] = faceZID_; 
    fieldType[JxID_] = edgeXID_; 
    fieldType[JyID_] = edgeYID_; 
    fieldType[JzID_] = edgeZID_; 
    fieldType[Bx_tm1ID_] = faceXID_; 
    fieldType[By_tm1ID_] = faceYID_; 
    fieldType[Bz_tm1ID_] = faceZID_; 
    fieldType[rhoID_] = vertID_; 

    return fieldType; 
};

/// deletes fieldType_ array 
void Grid::deleteFieldType_() { 
    delete [] fieldType_; 
};

/// constructs and returns fieldPtr__ array 
/*! fieldPtr_ is an nFieldsTotal_ array storing each field, 
 * so that they can be accessed via fieldID \n 
 * e.g. int fieldID = ExID_; \n 
 * double*** field = fieldPtr_[fieldID]; 
 */ 
double**** Grid::setFieldPtr_() { 
    double**** fieldPtr = new double***[nFieldsTotal_]; 
    
    fieldPtr[ExID_] = Ex_; 
    fieldPtr[EyID_] = Ey_; 
    fieldPtr[EzID_] = Ez_; 
    fieldPtr[BxID_] = Bx_; 
    fieldPtr[ByID_] = By_; 
    fieldPtr[BzID_] = Bz_; 
    fieldPtr[JxID_] = Jx_; 
    fieldPtr[JyID_] = Jy_; 
    fieldPtr[JzID_] = Jz_; 
    fieldPtr[Bx_tm1ID_] = Bx_tm1_; 
    fieldPtr[By_tm1ID_] = By_tm1_; 
    fieldPtr[Bz_tm1ID_] = Bz_tm1_; 
    fieldPtr[rhoID_] = rho_; 
   
    return fieldPtr; 
}; 

/// deletes fieldPtr_ array
void Grid::deleteFieldPtr_() { 
    delete [] fieldPtr_; 
}; 

/// checks validity of input parameters for Grid constructor 
/*! asserts necessary conditions on each input (mainly positivity of many parameters). Terminates program if inputs are incorrect.
 */ 
void Grid::checkInput_() { 
    assert(nx_ > 2*nGhosts_); // to guarantee there is at least 1 physical cell
    assert(ny_ > 2*nGhosts_); 
    assert(nz_ > 2*nGhosts_); 
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
    assert(dx_ > 0); 
    assert(dy_ > 0); 
    assert(dz_ > 0);
    assert(maxPointsInPlane_ > 0); 
    assert(nFieldsToSend_ == 9); 
    assert(ghostVecSize_ > 0); 
}; 

/// sets field corresponding to fieldID to zero
void Grid::zeroField_(const int fieldID) { 
    assert (fieldID > -1 && fieldID < nFieldsTotal_); 
    double*** field = fieldPtr_[fieldID]; 
    int i,j,k; // iterators 
    for (i=0; i<nxTot_; ++i) { 
        for (j=0; j<nyTot_; ++j) { 
            for (k=0; k<nzTot_; ++k) { 
                field[i][j][k]=0; 
            } 
        } 
    } 
} 

/// sets all components of J to be identically zero
void Grid::zeroJ() { 
    zeroField_(JxID_); 
    zeroField_(JyID_); 
    zeroField_(JzID_); 
};

/// sets rho to be identically zero
void Grid::zeroRho() { 
    zeroField_(rhoID_); 
};

/// sets all components of E to be identically zero
void Grid::zeroE() { 
    zeroField_(ExID_); 
    zeroField_(EyID_); 
    zeroField_(EzID_); 
};

/// sets all components of B and B_tm1 to be identically zero
void Grid::zeroB() { 
    zeroField_(BxID_); 
    zeroField_(ByID_); 
    zeroField_(BzID_); 
    zeroField_(Bx_tm1ID_); 
    zeroField_(By_tm1ID_); 
    zeroField_(Bz_tm1ID_); 
};

/// Initialize E and B fields
/*! Use restart file to set values of initial E,B,J fields
 */ 
void Grid::InitializeFields(void){

    if(rank_MPI==0)printf("        Initializing fields by reading files...\n");
    // placeholder until restart files exist 
    zeroB(); 
    zeroE(); 
    zeroJ(); 
    zeroRho(); 
}; 

//! Execute field boundary conditions
void Grid::executeBC(void){
    // loop through dimensions
    for(int i=0;i<3;i++){
        // left and right boundary in each dimension
        boundaries_[2*i]->completeBC();
        boundaries_[2*i+1]->completeBC();
    }
};
