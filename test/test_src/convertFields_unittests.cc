#define MAIN_CPP
#include "gtest/gtest.h"
#include "math.h"
#include "poisson.hpp"

// Tests of input/output to public methods
class ConvertTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    //MPI
    rank_MPI = 0;
    size_MPI = 1;
    //Read input file
    Input *input =  new Input();
    char filename[100];
    sprintf(filename, "../data/unitest/input.txt");
    input->readinfo(filename);
    Input_Info_t *input_info = input->getinfo();

    //Initialize a domain and a grid
    Domain *domain = new Domain(input_info->nCell, input_info->nProc, input_info->xyz0, input_info->Lxyz);
    grid = new Poisson_Solver(domain, input_info);
  }

  virtual void TearDown() {
    delete grid;
  }

  Poisson_Solver *grid;
};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


// tests requiring internal examination
class ConvertPrivateTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    //MPI
    rank_MPI = 0;
    size_MPI = 1;
    //Read input file
    Input *input =  new Input();
    char filename[100];
    sprintf(filename, "../data/unitest/input.txt");
    input->readinfo(filename);
    Input_Info_t *input_info = input->getinfo();

    //Initialize a domain and a grid
    Domain *domain = new Domain(input_info->nCell, input_info->nProc, input_info->xyz0, input_info->Lxyz);
    grid = new Poisson_Solver(domain, input_info);
  }

  virtual void TearDown() {
    delete grid;
  }

  double sumField(double*** field) { 
    int i,j,k;
    double fieldSum = 0; 
    for (i=0; i<grid->nxTot_; ++i) { 
      for (j=0; j<grid->nyTot_; ++j) { 
	for (k=0; k<grid->nzTot_; ++k) { 
	  fieldSum += field[i][j][k]; 
	}
      }
    }
    return fieldSum; 
  }; 

  Poisson_Solver *grid;
};

// test that conversion phi to E is correct
// case: constant phi 
TEST_F(ConvertPrivateTest, constantPhiTest) {

  // set phi1 to constant 
  int i,j,k; 
  double phival = 1.0; 
  for (i=0; i<grid->nxTot_; ++i) { 
    for (j=0; j<grid->nyTot_; ++j) { 
      for (k=0; k<grid->nzTot_; ++k) { 
	grid->phi1_[i][j][k] = phival; 
      } 
    } 
  } 

  // derive E
  grid->constE(0,0,0);
  grid->phiToE();

  // for constant phi, E should vanish everywhere 
  double  ExSum = sumField(grid->Ex_);
  EXPECT_EQ(ExSum,0);
}

