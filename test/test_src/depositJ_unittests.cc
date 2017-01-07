#define MAIN_CPP
#include <stdio.h>
#include "gtest/gtest.h"
#include "math.h"
#include "poisson.hpp"
#include "particle_handler.hpp"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// tests requiring internal examination
class DepositJTest : public ::testing::Test {
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
    input_info = input->getinfo();

    //Initialize the domain, grid, and particle_handler
    domain = new Domain(input_info->nCell, input_info->nProc, input_info->xyz0, input_info->Lxyz);
    grid = new Poisson_Solver(domain, input_info);
    part_handler = new Particle_Handler();
    part_handler->Load(input_info,domain);
  }

  virtual void TearDown() {
    delete grid;
  }

  Poisson_Solver *grid;
  Particle_Handler *part_handler;
  Domain *domain;
  Input_Info_t *input_info;
};

// Test particle current-deposition is working.
TEST_F(DepositJTest, sumOverJ) {

  std::vector<Particle> parts = part_handler->getParticleVector();

  double mod_v = pow(parts[0].v[0],2.0)+pow(parts[0].v[1],2.0)+pow(parts[0].v[2],2.0);
  printf("mod_v=%f\n",mod_v);

  part_handler->depositRhoJ(grid,false,domain,input_info);
  printf("Passed depositRhoJ");

  // sum over J
  int i,j,k;
  double JSum = 0.0;
  for (i=0; i<grid->nxTot_; ++i) {
    for (j=0; j<grid->nyTot_; ++j) {
      for (k=0; k<grid->nzTot_; ++k) {
        JSum += grid->Jx_[i][j][k];
	JSum += grid->Jy_[i][j][k];
	JSum += grid->Jz_[i][j][k];
      }
    }
  }

  //EXPECT_EQ(JSum,0);
  EXPECT_NEAR(JSum,0,0.1);
}

