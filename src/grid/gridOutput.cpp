#include "grid.hpp"
#include <assert.h>

/// get dimensions of physical region of field
void Grid::getDimPhys(const int fieldID, int* dim) {
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
    // directions 
    int xside=1; 
    int yside=2; 
    int zside=3; 
    
    // limits 
    int iEnd = sideToIndex_(xside,fieldID)+1;  
    int jEnd = sideToIndex_(yside,fieldID)+1; 
    int kEnd = sideToIndex_(zside,fieldID)+1; 
    

    // local dims for this field
    dim[0] = iEnd-iBeg_;
    dim[1] = jEnd-jBeg_;
    dim[2] = kEnd-kBeg_;
}

///// write all physical points of single field component to output 
//void Grid::writeFieldTimeseries_(FieldTimeseriesIO* tsIO, const int fieldID, const int iwrite) { 
//    // check for legal fieldID and make sure this field ID has an allocated hdf5 dataset 
//    assert(fieldID > -1 && fieldID < tsIO->getnFieldDatasets()); 
//
//    // use fieldID to get the pointer to the field 
//    double*** mat = fieldPtr_[fieldID]; 
//
//    tsIO->writeField(fieldID, mat, iwrite);
//} 
//
///// writes all physical points of J to output 
//void Grid::J_writeFieldTimeseries(FieldTimeseriesIO *tsio, const int iwrite) { 
//    writeFieldTimeseries_(tsio, JxID_, iwrite); 
//    writeFieldTimeseries_(tsio, JyID_, iwrite); 
//    writeFieldTimeseries_(tsio, JzID_, iwrite); 
//} 
//
///// writes all physical points of E to output 
//void Grid::E_writeFieldTimeseries(FieldTimeseriesIO *tsio, const int iwrite) { 
//    writeFieldTimeseries_(tsio, ExID_, iwrite); 
//    writeFieldTimeseries_(tsio, EyID_, iwrite); 
//    writeFieldTimeseries_(tsio, EzID_, iwrite);
//} 
//
///// writes all physical points of B to output 
//void Grid::B_writeFieldTimeseries(FieldTimeseriesIO *tsio, const int iwrite) { 
//    writeFieldTimeseries_(tsio, BxID_, iwrite); 
//    writeFieldTimeseries_(tsio, ByID_, iwrite); 
//    writeFieldTimeseries_(tsio, BzID_, iwrite); 
//} 
//
///// writes all physical points of rho to output 
//void Grid::Rho_writeFieldTimeseries(FieldTimeseriesIO *tsio, const int iwrite) { 
//    writeFieldTimeseries_(tsio, rhoID_, iwrite); 
//} 

/// slices a physical plane in the specified direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be written out. \n 
 * side is an integer -/+ 1 to indicate the location on the left/right side in the x direction, -/+ 2 in y, -/+ 3 in z. \n 
 * offset is an integer offset from the first/last physical index determined by side (e.g. side=-1 and offset=0 gives the yz plane of the 1st physical grid points in x direction, whereas offset=-1 would have returned the adjacent ghost cells and offset = 3 would have returned the 4th physical yz plane from the left). 
 */ 
<<<<<<< HEAD
void Grid::physSliceOut_(const int fieldID, const int side, const int offset) { 
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
    assert(side != 0 && abs(side) < 4); 

    // get the index to slice from
    int dex = sideToIndex_(side,fieldID) + offset; 
    assert(dex > -1); 
    
    // use fieldID to get the pointer to the field 
    double*** mat = fieldPtr_[fieldID]; 
   
    // directions 
    int xside=1; 
    int yside=2; 
    int zside=3; 
    
    // limits 
    int iEnd = sideToIndex_(xside,fieldID)+1;  
    int jEnd = sideToIndex_(yside,fieldID)+1; 
    int kEnd = sideToIndex_(zside,fieldID)+1; 
    
    // iterators 
    int i,j,k; 
    
    // slice an yz plane 
    if (abs(side)==1) {
        assert(dex < nxTot_);
        for (j=jBeg_; j<jEnd; ++j) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                // do stuff with hdf5 here with 
                // mat[dex][j][k]; 
                if (j==jBeg_ && k==kBeg_) {
                    printf("in physSliceOut_"); 
                } 
            }
        } 
    } 
    // slice an xz plane 
    else if (abs(side)==2) { 
        assert(dex < nyTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                // do exact same stuff with hdf5 here
                // note these cases are necessary because I'm not sure of any way 
                // to "loop" over these cases when they place dex in different 
                // orders of the [][][] notation 
                // mat[i][dex][k]; 
                if (i==iBeg_ && k==kBeg_) {
                    printf("in physSliceOut_"); 
                } 
            }
        } 
    } 
    // slice an xy plane 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (j=jBeg_; j<jEnd; ++j) { 
                // do exact same stuff with hdf5 here 
                // mat[i][j][dex]; 
                if (i==iBeg_ && j==jBeg_) {
                    printf("in physSliceOut_"); 
                } 
            }
        } 
    } 
}

/// writes a slice of physical points of J to output 
void Grid::JSliceOut(const int side, const int offset) { 
    physSliceOut_(JxID_,side,offset); 
    physSliceOut_(JyID_,side,offset); 
    physSliceOut_(JzID_,side,offset); 
} 

/// writes a slice of physical points of E to output 
void Grid::ESliceOut(const int side, const int offset) { 
    physSliceOut_(ExID_,side,offset); 
    physSliceOut_(EyID_,side,offset); 
    physSliceOut_(EzID_,side,offset);
} 

/// writes a slice of physical points of B to output 
void Grid::BSliceOut(const int side, const int offset) { 
    physSliceOut_(BxID_,side,offset); 
    physSliceOut_(ByID_,side,offset); 
    physSliceOut_(BzID_,side,offset); 
} 

/// writes a slice of  physical points of rho to output 
void Grid::RhoSliceOut(const int side, const int offset) { 
    physSliceOut_(rhoID_,side,offset); 
} 
=======
//void Grid::physSliceOut_(const int fieldID, const int side, const int offset) { 
//    // check for legal fieldID and side parameters 
//    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
//    assert(side != 0 && abs(side) < 4); 
//
//    // get the index to slice from
//    int dex = sideToIndex_(side,fieldID) + offset; 
//    assert(dex > -1); 
//    
//    // use fieldID to get the pointer to the field 
//    double*** mat = fieldPtr_[fieldID]; 
//   
//    // directions 
//    int xside=1; 
//    int yside=2; 
//    int zside=3; 
//    
//    // limits 
//    int iEnd = sideToIndex_(xside,fieldID)+1;  
//    int jEnd = sideToIndex_(yside,fieldID)+1; 
//    int kEnd = sideToIndex_(zside,fieldID)+1; 
//    
//    // iterators 
//    int i,j,k; 
//    
//    // slice an yz plane 
//    if (abs(side)==1) {
//        assert(dex < nxTot_);
//        for (j=jBeg_; j<jEnd; ++j) { 
//            for (k=kBeg_; k<kEnd; ++k) { 
//                // do stuff with hdf5 here
//                mat[dex][j][k]; 
//            }
//        } 
//    } 
//    // slice an xz plane 
//    else if (abs(side)==2) { 
//        assert(dex < nyTot_); 
//        for (i=iBeg_; i<iEnd; ++i) { 
//            for (k=kBeg_; k<kEnd; ++k) { 
//                // do exact same stuff with hdf5 here
//                // note these cases are necessary because I'm not sure of any way 
//                // to "loop" over these cases when they place dex in different 
//                // orders of the [][][] notation 
//                mat[i][dex][k]; 
//            }
//        } 
//    } 
//    // slice an xy plane 
//    else if (abs(side)==3) { 
//        assert(dex < nzTot_); 
//        for (i=iBeg_; i<iEnd; ++i) { 
//            for (j=jBeg_; j<jEnd; ++j) { 
//                // do exact same stuff with hdf5 here 
//                mat[i][j][dex]; 
//            }
//        } 
//    } 
//}
//
///// writes a slice of physical points of J to output 
//void Grid::JSliceOut(const int side, const int offset) { 
//    physSliceOut_(JxID_); 
//    physSliceOut_(JyID_); 
//    physSliceOut_(JzID_); 
//} 
//
///// writes a slice of physical points of E to output 
//void Grid::ESliceOut(const int side, const int offset) { 
//    physSliceOut_(ExID_); 
//    physSliceOut_(EyID_); 
//    physSliceOut_(EzID_);
//} 
//
///// writes a slice of physical points of B to output 
//void Grid::BSliceOut(const int side, const int offset) { 
//    physSliceOut_(BxID_); 
//    physSliceOut_(ByID_); 
//    physSliceOut_(BzID_); 
//} 
//
///// writes a slice of  physical points of rho to output 
//void Grid::RhoSliceOut(const int side, const int offset) { 
//    physSliceOut_(rhoID_); 
//} 

