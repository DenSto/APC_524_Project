#ifndef GRID_HPP
#define GRID_HPP

class Grid {
 public:
 	Grid(int nx, int ny, int nz, double x0, double y0, double z0, double dx);
	virtual ~Grid();

	int evolveFields (double dt);
    void updateGhostCells(); 
	int getFieldInterpolatorVec (int cellID, * double InterpolatorVec);
	int getCellID(double x, double y, double z);

	int getGhostVecSize(); // called by main to size MPI Buffer
	void getGhostVec(const int side, double* ghostVec); // called by main to get MPI 
	void setGhostVec(const int side, const double* ghostVec);

 protected:

 	const int nx_;     // number of (physical + ghost) gridpoints  
 	const int ny_;
 	const int nz_;

    const int nGhosts_; // number of ghost points in each dimension/2

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

    // vector for storing temporary physical slices of scalar fields
    double *sliceTmp_;

    void sliceMatToVec(const double*** mat, const int side, double* slice);
    void unsliceMatToVec(const double*** mat, const int side, double* lsice); 

 };

#endif
