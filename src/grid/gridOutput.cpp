/// write all physical points of single field component to output 
void Grid::physFieldOut_(const int fieldID) { 
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 

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
    for (i=iBeg_; i<iEnd; ++i) { 
        for (j=jBeg_; j<jEnd; ++j) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                
                // physical value at grid point i,j,k
                // probably do hdf5 stuff here 
                mat[i][j][k]; 
            }
        } 
    }
} 

/// writes all physical points of J to output 
void Grid::JFieldOut() { 
    physFieldOut_(JxID_); 
    physFieldOut_(JyID_); 
    physFieldOut_(JzID_); 
} 

/// writes all physical points of E to output 
void Grid::EFieldOut() { 
    physFieldOut_(ExID_); 
    physFieldOut_(EyID_); 
    physFieldOut_(EzID_);
} 

/// writes all physical points of B to output 
void Grid::BFieldOut() { 
    physFieldOut_(BxID_); 
    physFieldOut_(ByID_); 
    physFieldOut_(BzID_); 
} 

/// writes all physical points of rho to output 
void Grid::RhoFieldOut() { 
    physFieldOut_(rhoID_); 
} 

/// slices a physical plane in the specified direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be written out. \n 
 * side is an integer -/+ 1 to indicate the location on the left/right side in the x direction, -/+ 2 in y, -/+ 3 in z. \n 
 * offset is an integer offset from the first/last physical index determined by side (e.g. side=-1 and offset=0 gives the yz plane of the 1st physical grid points in x direction, whereas offset=-1 would have returned the adjacent ghost cells and offset = 3 would have returned the 4th physical yz plane from the left). 
 */ 
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
                // do stuff with hdf5 here
                mat[dex][j][k]; 
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
                mat[i][dex][k]; 
            }
        } 
    } 
    // slice an xy plane 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (j=jBeg_; j<jEnd; ++j) { 
                // do exact same stuff with hdf5 here 
                mat[i][j][dex]; 
            }
        } 
    } 
}

/// writes a slice of physical points of J to output 
void Grid::JSliceOut(const int side, const int offset) { 
    physSliceOut_(JxID_); 
    physSliceOut_(JyID_); 
    physSliceOut_(JzID_); 
} 

/// writes a slice of physical points of E to output 
void Grid::ESliceOut(const int side, const int offset) { 
    physSliceOut_(ExID_); 
    physSliceOut_(EyID_); 
    physSliceOut_(EzID_);
} 

/// writes a slice of physical points of B to output 
void Grid::BSliceOut(const int side, const int offset) { 
    physSliceOut_(BxID_); 
    physSliceOut_(ByID_); 
    physSliceOut_(BzID_); 
} 

/// writes a slice of  physical points of rho to output 
void Grid::RhoSliceOut(const int side, const int offset) { 
    physSliceOut_(rhoID_); 
} 

