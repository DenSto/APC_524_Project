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
#include <gtest/gtest_prod.h>
#include "../boundaries/fields_boundary.hpp"

class Grid {

public:
  Grid(int *nxyz, int nGhosts, double *xyz0, double *Lxyz);
  virtual ~Grid();

  int evolveFields (double dt);
  int evolveFieldsES (double dt);

  // Initialize fields by either solve Poisson's equation or read restart file.  
  virtual void InitializeFields(void);
  
  void constJ(double vx, double vy, double vz); 
  void constE(double vx, double vy, double vz); 
  void constB(double vx, double vy, double vz); 
  void constRho(double v); 
  
  int addJ(int cellID, double *Jvec);
  int addRho(int cellID, double *Rhovec);
  int getFieldInterpolatorVec (int cellID, double* InterpolatorVec);
  int getCellID(double x, double y, double z);
  int getCellVertex(int cellID, double *xyz);
  int getNumberOfCells();
  int getNumCells3D(double *nvec);
  double getStepSize(int dimension);

  int setFieldAlongEdge( std::string &fieldStr, int dim, bool edge, double fieldVal);

  virtual int getGhostVecSize(const int sendID); // called by main to size MPI Buffer
  // side = +1: x right, side = -1: x left
  // side = +2: y right, side = -2: y left
  // side = +3: z right, side = -3: z left
  virtual void getGhostVec(const int side, double* ghostVec, int sendID); // called by main to get MPI
  virtual void setGhostVec(const int side, double* ghostVec, int sendID);
  void updatePeriodicGhostCells(); 
  void setBoundaryVec(const int side, const double* ghostVec); // load physical boundary conditions
                                                               // boundary condition may depend on time_phys 

  void EFieldOut(); 
  void BFieldOut(); 
  void JFieldOut(); 
  void RhoFieldOut(); 
  
  void ESliceOut(const int side, const int offset); 
  void BSliceOut(const int side, const int offset); 
  void JSliceOut(const int side, const int offset); 
  void RhoSliceOut(const int side, const int offset); 
  
  void executeBC(void); // execute field boundary conditions
  void setBoundaries(BC_Field** bc){boundaries_=bc;}
  void freeBoundaries(void){delete [] boundaries_;}

  // need to be public, for use by field boundary conditions 
  double* sliceTmp_; 
  double* ghostTmp_;  

protected:
  BC_Field** boundaries_; // field boundary conditions

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

  const double dx_;
  const double dy_;
  const double dz_;

  const double idx_;
  const double idy_;
  const double idz_;

  const int maxPointsInPlane_;
  
  const int nFieldsTotal_;  
  const int nFieldsJEB_;
  const int ExID_; 
  const int EyID_; 
  const int EzID_; 
  const int BxID_; 
  const int ByID_; 
  const int BzID_; 
  const int JxID_; 
  const int JyID_; 
  const int JzID_; 
  const int Bx_tm1ID_; 
  const int By_tm1ID_; 
  const int Bz_tm1ID_;
  const int rhoID_; 
  
  const int nTypes_; 
  const int edgeXID_; 
  const int edgeYID_; 
  const int edgeZID_; 
  const int faceXID_; 
  const int faceYID_; 
  const int faceZID_; 
  const int vertID_; 

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
  double ***rho_; 

  int *fieldType_; 
  int **fieldSize_;  
  double ****fieldPtr_; 

  int *fieldIsContiguous_;
  

  // allocates contiguous memory for nx*ny*nz array
  double*** newField_(int ifield);
  // deallocates memory for nx*ny*nz array
  void deleteField_(double*** fieldPt,int ifield);
  int** setFieldSize_(); 
  void deleteFieldSize_(); 
  int* setFieldType_(); 
  void deleteFieldType_();
  double**** setFieldPtr_(); 
  void deleteFieldPtr_(); 

  void constField_(const int fieldID, const double val); 

  int sideToIndex_(const int side, const int fieldID);
  /* assert statements to check necessary conditions for initialized variables */
  void checkInput_();

  // stores a 2D plane of ghost points in sliceTmp_
  void sliceMatToVec_(const int fieldID, const int side, const int offset, double* vec);
  // puts a 2D plane of ghost points from sliceTmp_ into mat
  void unsliceMatToVec_(const int fieldID, const int side, const int offset, double* vec);

  int setFieldInPlane_( int dim, int indx, double *** field, double fieldVal);

  void physFieldOut_(const int fieldID);
  void physSliceOut_(const int fieldID, const int side, const int offset); 


  // for unit testing
  friend class oGridInternalTest;
  FRIEND_TEST( oGridInternalTest, EMWave);
  FRIEND_TEST( oGridInternalTest, EMWaveLong);

  // unit testing in grid_unittests.cc 
  friend class GridPrivateTest; 
  FRIEND_TEST(GridPrivateTest, fieldSizeTest); 
  FRIEND_TEST(GridPrivateTest, fieldPtrTest); 
  FRIEND_TEST(GridPrivateTest, zeroFields); 

  // unit testing in spookyGrid_unittests.cc
  FRIEND_TEST(GridPrivateTest, sideToIndexTest); 
  FRIEND_TEST(GridPrivateTest, periodicUpdateTest); 
  FRIEND_TEST(GridPrivateTest, ghostVecSizeTest); 

 };
#endif
