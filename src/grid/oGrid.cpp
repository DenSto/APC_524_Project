#include <cstdio>
#include <algorithm>
#include <math.h>
#include "grid.hpp"



//! Evolve Electric and Magnetic fields in time
/*!
	Uses Yee algorithm to advance E and B fields.
	Assumes Gaussian-style Maxwell equation, with c = 1. 
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
								  + 4*M_PI*Jx_[ix][iy][iz]  );
				Ey_[ix][iy][iz] = Ey_[ix][iy][iz] - dt * 
								  (  ( Bz_tm1_[ix][iy][iz] - Bz_tm1_[ix-1][iy][iz] )/dx_  
								   - ( Bx_tm1_[ix][iy][iz] - Bx_tm1_[ix][iy][iz-1] )/dz_ 
								 + 4*M_PI*Jy_[ix][iy][iz]  );
				Ez_[ix][iy][iz] = Ez_[ix][iy][iz] - dt * 
								  (  ( Bx_tm1_[ix][iy][iz] - Bx_tm1_[ix][iy-1][iz] )/dy_  
								   - ( By_tm1_[ix][iy][iz] - By_tm1_[ix-1][iy][iz] )/dx_ 
								  + 4*M_PI*Jz_[ix][iy][iz]  );
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


//! Add currents from particle to grid
/*!
	Currents added to cell with ID cellID via input vector of form:\n
	[Jx((0,0,0) -> (1,0,0)), Jx((0,1,0) -> (1,1,0)), Jx((0,1,1) -> (1,1,1)), Jx((0,0,1) -> (1,0,1)),...\n
	Jy((0,0,0) -> (0,1,0)), Jy((0,0,1) -> (0,1,1)), Jy((1,0,1) -> (1,1,1)), Jy((1,0,0) -> (1,1,0)),...\n
	Jz((0,0,0) -> (0,0,1)), Jz((1,0,0) -> (1,0,1)), Jz((1,1,0) -> (1,1,1)), Jz((0,1,0) -> (0,1,1))]
*/
int Grid::addJ(int cellID, double *Jvec) {
	// get indices for cellID
	int iz = cellID % nz_;
	int iy = ( cellID - iz) % (ny_*nz_);
	int ix = ( cellID - iz - nz_*iy ) / (ny_ * nz_);

	// put down currents
	Jx_[ix][iy][iz] += Jvec[0];
	Jx_[ix][iy+1][iz] += Jvec[1];
	Jx_[ix][iy+1][iz+1] += Jvec[2];
	Jx_[ix][iy+1][iz] += Jvec[3];

	Jy_[ix][iy][iz] += Jvec[4];
	Jy_[ix][iy][iz+1] += Jvec[5];
	Jy_[ix+1][iy][iz+1] += Jvec[6];
	Jy_[ix+1][iy][iz] += Jvec[7];

	Jz_[ix][iy][iz] += Jvec[8];
	Jz_[ix+1][iy][iz] += Jvec[9];
	Jz_[ix+1][iy+1][iz] += Jvec[10];
	Jz_[ix][iy+1][iz] += Jvec[11];

	return 0;
};



//! Return vector for field interpolation
/*!
	Based on cellID, return relevant edge E and face B fields and cell origin, in format:\n
	[x, y, z, ...\n
	Ex( ix, iy, iz ), Ex( ix, iy+1,iz ), Ex( ix, iy+1, iz+1 ), Ex( ix, iy, iz+1 ), ...\n
	Ey( ix, iy, iz ), Ey( ix, iy, iz+1 ), Ey( ix+1, iy, iz+1 ), Ey( ix+1, iy, iz ), ...\n
	Ez( ix, iy, iz ), Ez( ix+1, iy, iz ), Ez( ix+1, iy+1, iz ), Ez( ix, iy+1, iz ), ...\n
	Bx( ix, iy, iz ), Bx( ix+1, iy, iz ), ...\n
	By( ix, iy, iz ), By( ix, iy+1, iz ), ...\n
	Bz( ix, iy, iz ), Bz( ix, iy, iz+1 ), ...]\n
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
	\n
	If particle is in a ghost cell or off the grid entirely, returns \n
	-1 if off (-z), -2 if off (+z) \n
	-3 if off (-y), -4 if off (+y) \n
	-5 if off (-x), -6 if off (+x)
*/
int Grid::getCellID(double x, double y, double z) {
	// get indices in x, y, z
	int ix = ((int) (x-x0_))*idx_;
	int iy = ((int) (y-y0_))*idy_;
	int iz = ((int) (z-z0_))*idz_;

	// inform if in ghost cell, and which one
	if ( iz < nGhosts_ ) {
		return -1;
	} else if ( iz > nz_ - nGhosts_ - 1 ) {
		return -2;
	} else if ( iy < nGhosts_ ) {
		return -3; 
	} else if ( iy > ny_ - nGhosts_ - 1) {
		return -4;
	} else if ( ix < nGhosts_ ) {
		return -5;
	} else if ( ix > nx_ - nGhosts_ - 1) {
		return -6;
	}

	return (ny_*nz_)*ix + nz_*iy + iz;

};

//! Returns vertex corresponding to cell ID
/*
	Vertex stored in input vector *xyz.
*/
int Grid::getCellVertex(int cellID, double *xyz) {
	int iz = cellID % nz_;
	int iy = ( cellID - iz) % (ny_*nz_);
	int ix = ( cellID - iz - nz_*iy ) / (ny_ * nz_);

	// x, y, z
	xyz[0] = x0_ + dx_ * ix;
	xyz[1] = y0_ + dy_ * iy;
	xyz[2] = z0_ + dz_ * iz;

	return 0;
};

//! Get total number of cells in grid.
/*!
	Includes ghost cells.
*/
int Grid::getNumberOfCells() {
	return nx_*ny_*nz_;
};

//! Get step size along dimension in grid.
/*!
	Returns step size along dimension according to;
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
		default: printf("Invalid dimension in Grid::getStepSize\n"); return -1;
	}
};

