#define MAIN_CPP
#include "gtest/gtest.h"
#include "domain.hpp"
#include "grid.hpp"
#include "hdf5io.hpp"
#include "math.h"

// Tests of input/output to public methods
class FieldIOTest : public ::testing::Test {
protected:
  virtual void SetUp() {

  }

  virtual void TearDown() {
    delete grid;
  }

  Domain *domain;
  Grid *grid;
  Hdf5IO *hdf5io;
  FieldTimeseriesIO *field_tsio;
};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST_F(FieldIOTest, writeField1D) {
    int nxyz [3] = {3,3,4};
    int nProc [3] = {1,1,1};
    int nGhosts = 1; 
    double xyz0 [3] = {0,0,0};
    double Lxyz [3] = {1,1,1};

    int which_fields = 0; //only write rho
    int nwrite = 2;

    size_MPI = 1;
    rank_MPI = 0;
    domain = new Domain(nxyz, nProc, xyz0, Lxyz);
    grid = new Grid(nxyz, nGhosts, xyz0, Lxyz);
    hdf5io = new Hdf5IO("test_data/hdf5io_unittests.h5");
    field_tsio = new FieldTimeseriesIO(hdf5io, grid, domain, which_fields, nwrite);

    grid->constRho(0.);
    field_tsio->writeFields(grid, which_fields, 0);

    int i,j,k; 
    for (i=0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                grid->rho_[i][j][k] = k + 10*j + 100*i; 
		printf("rho(%d,%d,%d) = %f\n", i,j,k, grid->rho_[i][j][k]);
            }
        }
    }
    field_tsio->writeFields(grid, which_fields, 1);


}


//
//// tests requiring internal examination
//class GridPrivateTest : public ::testing::Test {
//protected:
//  virtual void SetUp() {
//    int nxyz [3] = {3,5,7};
//    int nGhosts = 1; 
//    double xyz0 [3] = {0,0,0};
//    double Lxyz [3] = {1,1,1};
//
//    grid = new Grid(nxyz, nGhosts, xyz0, Lxyz);
//  }
//
//  virtual void TearDown() {
//    delete grid;
//  }
//
//  double sumField(double*** field) { 
//    int i,j,k;
//    double fieldSum = 0; 
//    for (i=0; i<grid->nxTot_; ++i) { 
//      for (j=0; j<grid->nyTot_; ++j) { 
//	for (k=0; k<grid->nzTot_; ++k) { 
//	  fieldSum += field[i][j][k]; 
//	}
//      }
//    }
//    return fieldSum; 
//  }; 
//
//  Grid *grid;
//
//};
//
//// test that fieldSize is correct  
//TEST_F(GridPrivateTest, fieldSizeTest) {
//    int xdir=0; 
//    int ydir=1; 
//    int zdir=2; 
//
//    // test a set of edges
//    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][xdir],grid->nx_); 
//    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][ydir],grid->nyTot_); 
//    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][zdir],grid->nzTot_); 
//
//    // test a set of faces 
//    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][xdir],grid->nx_); 
//    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][ydir],grid->nyTot_); 
//    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][zdir],grid->nz_);
//    
//    // test a set of vertices
//    EXPECT_EQ(grid->fieldSize_[grid->vertID_][xdir],grid->nxTot_); 
//    EXPECT_EQ(grid->fieldSize_[grid->vertID_][ydir],grid->nyTot_); 
//    EXPECT_EQ(grid->fieldSize_[grid->vertID_][zdir],grid->nzTot_); 
//}
//
//// test that fieldPtr works 
//TEST_F(GridPrivateTest, fieldPtrTest) {
//    double Exval = 1; 
//    double Byval = 2; 
//    double rhoval = 3; 
//
//    int i,j,k; 
//    for (i=0; i<grid->nxTot_; ++i) { 
//        for (j=0; j<grid->nyTot_; ++j) { 
//            for (k=0; k<grid->nzTot_; ++k) { 
//                grid->Ex_[i][j][k] = Exval; 
//                grid->By_[i][j][k] = Byval; 
//                grid->rho_[i][j][k] = rhoval; 
//            }
//        }
//    }
//
//    double*** ExPtr = grid->fieldPtr_[grid->ExID_]; 
//    double*** ByPtr = grid->fieldPtr_[grid->ByID_];
//    double*** rhoPtr = grid->fieldPtr_[grid->rhoID_];
//
//    EXPECT_EQ(grid->Ex_[1][1][1],ExPtr[1][1][1]); 
//    EXPECT_EQ(grid->By_[1][1][1],ByPtr[1][1][1]); 
//    EXPECT_EQ(grid->rho_[1][1][1],rhoPtr[1][1][1]); 
//}
//
//
//// test that zeroField methods work 
//TEST_F(GridPrivateTest,zeroFields) {
// 
//    // zero all fields 
//    grid->constE(0,0,0); 
//    grid->constB(0,0,0); 
//    grid->constJ(0,0,0); 
//    grid->constRho(0); 
//
//    // sum each field 
//    double ExSum = sumField(grid->Ex_); 
//    double BySum = sumField(grid->By_); 
//    double JzSum = sumField(grid->Jz_); 
//    double rhoSum = sumField(grid->rho_); 
//
//    // verify that each field sums to 0
//    EXPECT_EQ(ExSum,0); 
//    EXPECT_EQ(BySum,0); 
//    EXPECT_EQ(JzSum,0); 
//    EXPECT_EQ(rhoSum,0); 
//
//    // set one point in Ex to this value 
//    double addval = -2.0; 
//    grid->Ex_[1][1][1] += addval; 
//    ExSum = sumField(grid->Ex_); 
//    // verify that the sum is no longer zero 
//    EXPECT_EQ(ExSum,addval); 
//    // zero E again 
//    grid->constE(0,0,0); 
//    // verify that E again sums to zero 
//    ExSum = sumField(grid->Ex_); 
//    EXPECT_EQ(ExSum,0); 
//
//}




