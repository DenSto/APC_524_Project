#include <cstdio>
#include <algorithm>
#include "grid.hpp"

//! Evolve Electric and Magnetic fields in time
/*!
	Uses Yee algorithm to advance E and B fields.
*/
int Grid::evolveFields (double dt) {

	// set B(t-1) = B(t)
	std::swap(Bx_tm1_, Bx_);
	std::swap(By_tm1_, By_);
	std::swap(Bz_tm1_, Bz_);

	// calculate E 
	for (int ix = 1; ix < nx_; ix++) {
		for (int iy = 1; iy < ny_; iy++) {
			for (int iz = 1; iz < nz_; iz++) {
				Ex_[ix][iy][iz] = Ex_[ix][iy][iz] - dt * 
								  (  ( By_tm1_[ix][iy][iz] - By_tm1_[ix][iy][iz-1] )/dz_  
								   - ( Bz_tm1_[ix][iy][iz] - Bz_tm1_[ix][iy-1][iz] )/dy_ 
								  + Jx_[ix][iy][iz]  );
				Ey_[ix][iy][iz] = Ey_[ix][iy][iz] - dt * 
								  (  ( Bz_tm1_[ix][iy][iz] - Bz_tm1_[ix-1][iy][iz] )/dx_  
								   - ( Bx_tm1_[ix][iy][iz] - Bx_tm1_[ix][iy][iz-1] )/dz_ 
								 + Jy_[ix][iy][iz]  );
				Ez_[ix][iy][iz] = Ez_[ix][iy][iz] - dt * 
								  (  ( Bx_tm1_[ix][iy][iz] - Bx_tm1_[ix][iy-1][iz] )/dy_  
								   - ( By_tm1_[ix][iy][iz] - By_tm1_[ix-1][iy][iz] )/dx_ 
								  + Jz_[ix][iy][iz]  );
			}
		}
	}

	// calculate B
	for (int ix = 1; ix < nx_; ix++) {
		for (int iy = 1; iy < ny_; iy++) {
			for (int iz = 1; iz < nz_; iz++) {
				Bx_[ix][iy][iz] = Bx_[ix][iy][iz] + dt * 
								  (  ( Ey_[ix][iy][iz+1] - Ey_[ix][iy][iz] )/dz_  
								   - ( Ez_[ix][iy+1][iz] - Ez_[ix][iy][iz] )/dy_ );
				By_[ix][iy][iz] = By_[ix][iy][iz] + dt * 
								  (  ( Ez_[ix+1][iy][iz] - Ez_[ix][iy][iz] )/dx_  
								   - ( Ex_[ix][iy][iz+1] - Ex_[ix][iy][iz] )/dz_ );
				Bz_[ix][iy][iz] = Bz_[ix][iy][iz] + dt * 
								  (  ( Ex_[ix][iy+1][iz] - Ex_[ix][iy][iz] )/dy_  
								   - ( Ey_[ix+1][iy][iz] - Ey_[ix][iy][iz] )/dx_ );
			}
		}
	}

	return 0;
};

//! Return vector for field interpolation
/*!
	Based on cellID, return relevant edge E and face B fields and cell origin, in format
	[x, y, z, ...
	Ex( ix, iy, iz ), Ex( ix, iy+1,iz ), Ex( ix, iy+1, iz+1 ), Ex( ix, iy, iz+1 ), ...
	Ey( ix, iy, iz ), Ey( ix, iy, iz+1 ), Ey( ix+1, iy, iz+1 ), Ey( ix+1, iy, iz ), ...
	Ez( ix, iy, iz ), Ez( ix+1, iy, iz ), Ez( ix+1, iy+1, iz ), Ez( ix, iy+1, iz ), ...
	Bx( ix, iy, iz ), Bx( ix+1, iy, iz ), ...
	By( ix, iy, iz ), By( ix, iy+1, iz ), ...
	Bz( ix, iy, iz ), Bz( ix, iy, iz+1 ), ...]
	where ix, iy, and iz are the row indices for each of the three dimensions (calculated from the cellID)
*/
int Grid::getFieldInterpolatorVec (int cellID, double* InterpolatorVec) {
	//invert from cellID to get indices in x, y, z
	int iz = cellID % nz_;
	int iy = ( cellID - iz) % (ny_*nz_);
	int ix = ( cellID - iz - nz_*iy ) / (ny_ * nz_);

	// x, y, z
	InterpolatorVec[0] = x0_ + dx_ * ix;
	InterpolatorVec[1] = y0_ + dy_ * iy;
	InterpolatorVec[2] = z0_ + dz_ * iz;

	//Ex
	InterpolatorVec[3] = Ex_[ix][iy][iz];
	InterpolatorVec[4] = Ex_[ix][iy+1][iz];
	InterpolatorVec[5] = Ex_[ix][iy+1][iz+1];
	InterpolatorVec[6] = Ex_[ix][iy][iz+1];

	//Ey
	InterpolatorVec[7] = Ey_[ix][iy][iz];
	InterpolatorVec[8] = Ey_[ix][iy][iz+1];
	InterpolatorVec[9] = Ey_[ix+1][iy][iz+1];
	InterpolatorVec[10] = Ey_[ix+1][iy][iz];

	//Ez
	InterpolatorVec[11] = Ez_[ix][iy][iz];
	InterpolatorVec[12] = Ez_[ix+1][iy][iz];
	InterpolatorVec[13] = Ez_[ix+1][iy+1][iz];
	InterpolatorVec[14] = Ez_[ix][iy+1][iz];

	//Bx
	InterpolatorVec[15] = (Bx_tm1_[ix][iy][iz] + Bx_[ix][iy][iz])/2;
	InterpolatorVec[16] = (Bx_tm1_[ix+1][iy][iz] + Bx_[ix+1][iy][iz])/2;

	//By
	InterpolatorVec[17] = (By_tm1_[ix][iy][iz] + By_[ix][iy][iz])/2;
	InterpolatorVec[18] = (By_tm1_[ix][iy+1][iz] + By_[ix][iy+1][iz])/2;

	//Bz
	InterpolatorVec[19] = (Bz_tm1_[ix][iy][iz] + Bz_[ix][iy][iz])/2;
	InterpolatorVec[20] = (Bz_tm1_[ix][iy][iz+1] + Bz_[ix][iy][iz+1])/2;

	return 0;
};

//! Get cell ID based on particle position.
/*!
	Cell ID is uniquely given by (ny_*nz_)*ix + nz_*iy + iz.
	Returns -1 if particle is not on grid.
*/
int Grid::getCellID(double x, double y, double z) {
	//get indices in x, y, z
	int ix = ((int) (x-x0_))/dx_;
	int iy = ((int) (y-y0_))/dy_;
	int iz = ((int) (z-z0_))/dz_;

	if ( ( ix < 0 ) || ( ix > nx_-1 ) || ( iy < 0 ) || ( iy > ny_-1 ) || ( iz < 0 ) || ( iz > nz_-1 ) ) {
		printf("Particle out of grid bounds in Grid::getCellID\n");
		return -1;
	}

	return (ny_*nz_)*ix + nz_*iy + iz;

};

//! Get total number of cells.
int Grid::getNumberOfCells() {
	return nx_*ny_*nz_;
};

//! Get number of cells along dimension in grid.
/*!
	Returns number of cells along dimension according to;
	dimension = 0: x
	dimension = 1: y
	dimension = 2: z
	Returns -1 if invalid dimension.
*/
double Grid::getStepSize(int dimension) {
	switch(dimension) {
		case 0: return dx_; break;
		case 1: return dy_; break;
		case 2: return dz_; break;
		default: printf("Invalid dimension in Grid::getNumberOfCells\n"); return -1;
	}
};

