#define MAIN_CPP
#include <stdio.h>
#include "gtest/gtest.h"
#include "math.h"
#include "poisson.hpp"
#include "../boundaries/field_bc_factory.hpp"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// tests requiring internal examination
class PoissonTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    //MPI
    rank_MPI = 0;
    size_MPI = 1;
    //Read input file
    Input *input =  new Input();
    char filename[200];
    sprintf(filename, "test_data/poisson_unittest_data.txt");
    input->readinfo(filename);
    input_info = input->getinfo();

    //Initialize the domain and grid
    domain = new Domain(input_info->nCell, input_info->nProc, input_info->xyz0, input_info->Lxyz);
    grid = new Poisson_Solver(domain, input_info);

    //Initialize the field boundary conditions
    Field_BC_Factory::getInstance().Construct(domain,grid,input_info);
  }

  virtual void TearDown() {
    delete grid;
  }

  Poisson_Solver *grid;
  Domain *domain;
  Input_Info_t *input_info;
};

// Test field interpolation is working.
TEST_F(PoissonTest, testPoisson) {
  //Make const rho field rho = 0;
  grid->constRho(0.0);
  grid->constPhi(0.0);
  //Add test charges <--remember Poisson requires neutral charge distribution
  grid->rho_[3][3][3]=2.0; // cation at P=(3,3,3)
  grid->rho_[3][3][4]=-2.0;// anion at Q=(3,3,4)

  //Solve poisson's equation to high convergence.
  grid->run_poisson_solver_(grid->phi1ID_,grid->phi2ID_,grid->phi1_,grid->phi2_,grid->rho_,0.000001,-4*3.1415936535898);
  printf("Poisson converged!\n");

  //Get cell lengths
  double lcell[3] = {};
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  //Check point chage scales as 1/r  <--Only checks in z direction, which was made very long for this test to reduce wrap-around effects.
  double dz = lcell[2]; //dz
  double phiTest1 = grid->phi1_[3][3][10]; //phi(P+7*dz) or phi(Q+8*dz)
  double rP1 = 7.0*dz;
  double rQ1 = 8.0*dz;
  double phiTest2 = grid->phi1_[3][3][11]; //phi(P+8*dz) or phi(Q+9*dz)
  double rP2 = 8.0*dz;
  double rQ2 = 9.0*dz;

  EXPECT_NEAR(phiTest1/(1/rP1-1/rQ1),phiTest2/(1/rP2-1/rQ2),0.0001);
}
