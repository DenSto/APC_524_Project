#include <grid.hpp> 

Grid::Grid(int nx, int ny, int nz, double x0, double y0, double z0, double dx) { 
    
    nGhosts_ = 1; 
    sliceTmp_ = new double[nx_ - 2*nGhosts_]; 
} 

Grid::~Grid() { 
    
    delete [] sliceTmp_; 
} 

	int evolveFields (double dt);
	int getFieldInterpolatorVec (int cellID, * double InterpolatorVec);
	int getCellID(double x, double y, double z);

 	const int nx_;     // number of (physical + ghost) gridpoints  
 	const int ny_;
 	const int nz_;

    const int nGhosts_; // number of ghost points in each direction

    const double x0_;	// initial x position
 	const double y0_;	// initial y position
 	const double z0_;	// initial z position
    
 	const double dx_;

 	double ***Ex_;
 	double ***Ey_;
 	double ***Ez_;

 	double ***Bx_;
 	double ***By_;
 	double ***Bz_;

  	double ***Bx_tm1_; //timestep back (since B at half timesteps)
 	double ***By_tm1_;
 	double ***Bz_tm1_;

  	double ***Jx_;
 	double ***Jy_;
 	double ***Jz_;
