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
  //Add a test charge at grid pt (3,3,3)
  grid->rho_[3][3][3]=2.0;

  //Solve poisson's equation to high convergence.
  grid->run_poisson_solver_(grid->phi1ID_,grid->phi1_,grid->phi2_,grid->rho_,0.000000001,-4*3.1415936535898);

  //Get cell lengths
  double lcell[3] = {};
  for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

  //Check point chage scales as 1/r
  //double phi1 = grid->phi1_[3][3][3]; //at test charge
  double phi2 = grid->phi1_[3][3][4]; //dz away
  double phi3 = grid->phi1_[3][3][5]; //2*dz away
  double phi4 = grid->phi1_[2][3][3]; //dx away
  double r12 = lcell[2]; //dz
  double r13 = 2.0*lcell[2]; //2*dz
  double r14 = lcell[0]; //dx

  EXPECT_NEAR(phi2*r12,phi3*r13,0.1);
  EXPECT_NEAR(phi2*r12,phi4*r14,0.1);
}
