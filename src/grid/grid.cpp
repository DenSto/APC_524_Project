#include <grid.hpp> 

Grid::Grid(int nx, int ny, int nz, double x0, double y0, double z0, double dx) { 
    
    nGhosts_ = 1; 
    sliceTmp_ = new double[nx_ - 2*nGhosts_]; 
} 

Grid::~Grid() { 
    
    delete [] sliceTmp_; 
} 
