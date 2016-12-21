#ifndef GRID_HPP
#define GRID_HPP

//! Class representing grid on which E and B fields and currents are defined.
/*!
	Grid has ghost cells on each face.
	The ghost cell updating in y and z arises from periodic boundary conditions.
	x-direction ghost cells allow communication between MPI domains.

	Following Yee (1966), electric fields and currents reside on edges, and magnetic fields on faces.
	Fields are updated using a set of finite-difference equations approximating Ampere's and Faraday's Laws.

	A set of getters are available to allow particles to interpolate electric fields based on their position.

*/
#include <string>

class Grid {

public:
  Grid(int *nxyz, int nGhosts, double *xyz0, double *Lxyz);
  virtual ~Grid();

  int evolveFields (double dt);

  void InitializeFields(int restart); // Initialize E field by either solve Poisson's equation
                                      // or read restart file. Initialize B field with prescribed
                                      // geometry or read restart file. 
  void zeroJ();
  void zeroRho();
  int addJ(int cellID, double *Jvec);
  int addRho(int cellID, double *Rhovec);
  int getFieldInterpolatorVec (int cellID, double* InterpolatorVec);
  int getCellID(double x, double y, double z);
  int getCellVertex(int cellID, double *xyz);
  int getNumberOfCells();
  double getStepSize(int dimension);

  int setFieldAlongEdge( std::string &fieldStr, int dim, bool edge, double fieldVal);

  void updatePeriodicGhostCells();
  int getGhostVecSize(); // called by main to size MPI Buffer
  // side = +1: x right, side = -1: x left
  // side = +2: y right, side = -2: y left
  // side = +3: z right, side = -3: z left
  void getGhostVec(const int side, double* ghostVec); // called by main to get MPI
  void getGhostVecAlt(const int side, double* ghostVec); // called by main to get MPI
  void setGhostVec(const int side, const double* ghostVec);
  void setGhostVecAlt(const int side, const double* ghostVec);
  void setBoundaryVec(const int side, const double* ghostVec); // load physical boundary conditions
                                                               // boundary condition may depend on time_phys 

protected:

  const int nx_;     // number of (physical + ghost) gridpoints
  const int ny_;
  const int nz_;

  const int nGhosts_; // number of ghost points in each dimension/2

  const int nxTot_; // nx_ + 1 (total number of grid points in x) 
  const int nyTot_; 
  const int nzTot_; 

  const double x0_;	// initial x position
  const double y0_;	// initial y position
  const double z0_;	// initial z position

  const double Lx_;
  const double Ly_;
  const double Lz_;

  const int iBeg_; // indices marking beginning and end of physical
  const int jBeg_; // (non ghost) points in each direction
  const int kBeg_;
  const int iEnd_;
  const int jEnd_;
  const int kEnd_; 

  const double dx_;
  const double dy_;
  const double dz_;

  const double idx_;
  const double idy_;
  const double idz_;

  const int nRealPtsYZPlane_;
  const int nFieldsToSend_;
  const int nFieldsTotal_; 
  const int ghostVecSize_; /* total number of ghost field values in
                                    a single plane. All MPI communiation
                                    of fields send messages of this size */

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

  double ***rhox_; 
  double ***rhoy_; 
  double ***rhoz_; 

  // vector for storing temporary physical slices of scalar fields
  double *sliceTmp_;

  double *fieldIsContiguous_; 

  // allocates contiguous memory for nx*ny*nz array
  double*** newField_(int ifield);
  // deallocates memory for nx*ny*nz array
  void deleteField_(double*** fieldPt,int ifield);

  // converts side = -/+ 1 into a real index
  int sideToIndex_(const int side);
  /* assert statements to check necessary conditions for initialized variables */
  void checkInput_();

  // stores a 2D plane of ghost points in sliceTmp_
  void sliceMatToVec_(double*** const mat, const int side);
  // puts a 2D plane of ghost points from sliceTmp_ into mat
  void unsliceMatToVec_(double*** mat, const int side);

  int setFieldInPlane_( int dim, int indx, double *** field, double fieldVal);

 };
#endif
