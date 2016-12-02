#ifndef GRID_HPP
#define GRID_HPP

class Grid {
 public:
 	Grid(int nx, int ny, int nz, double x0, double y0, double z0, double dx);
	virtual ~Grid();

	int evolveFields (double dt);
	int getFieldInterpolatorVec (int cellID, * double InterpolatorVec);
	int getCellID(double x, double y, double z);

	int getGhostVecSize(); // called by main to size MPI Buffer
	int getGhostVec(int side, * double ghostVec); // called by main to get MPI 
	int setGhostVec(int side, const * double ghostVec);




 protected:
 	
 	const int nx_;     // number of gridpoints
 	const int ny_;
 	const int nz_;
 	const double x0_;	// initial x position
 	const double y0_;	// initial y position
 	const double z0_;	// initial y position
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


 };

#endif
