#include <grid.hpp> 

Grid::Grid(int nx, int ny, int nz, double x0, double y0, double z0, double dx) { 
    
    nx_ = nx; 
    ny_ = ny; 
    nz_ = nz; 

    nGhosts_ = 1; 
    ghostVecSize_ = 9*(ny_ - 2*nGhosts_)*(nz_ - 2*nGhosts_); 

    iBeg_ = 1; 
    jBeg_ = 1; 
    kBeg_ = 1; 
    iEnd_ = nx_ - (nGhosts_ + 1); 
    jEnd_ = ny_ - (nGhosts_ + 1); 
    kEnd_ = nz_ - (nGhosts_ + 1); 

    sliceTmp_ = new double[nx_ - 2*nGhosts_]; 
} 

Grid::~Grid() { 
    
    delete [] sliceTmp_; 
} 
