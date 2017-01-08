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

    int xgsize = grid->getGhostVecSize(-1);
    int ygsize = 1; //dummy
    int zgsize = 1; //dummy
    domain->mallocGhosts(xgsize,ygsize,zgsize);
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
  int np = part_handler->nParticles();
  double mod_v;
  double q;

  if (np != 1) {
    printf("This test did not fail, but can only be run for 1 particle. Change input file and rerun.");
  } else {
    //Get full vector of particles.
    std::vector<Particle> parts = part_handler->getParticleVector();

    //Calculate particle |v| and q
    mod_v = pow(pow(parts[0].v[0],2.0)+pow(parts[0].v[1],2.0)+pow(parts[0].v[2],2.0),0.5);
    q = parts[0].q;

    //Perform the deposition of the current (note: depositRho=false)
    part_handler->depositRhoJ(grid,false,domain,input_info);

    //Get cell volume
    double lcell[3] = {};
    for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);
    double cell_volume = lcell[0]*lcell[1]*lcell[2];

    // sum up all cells' J
    int i,j,k;
    double JxSum = 0.0;
    double JySum = 0.0;
    double JzSum = 0.0;
    for (i=0; i<grid->nxTot_; ++i) {
      for (j=0; j<grid->nyTot_; ++j) {
	for (k=0; k<grid->nzTot_; ++k) {
	  JxSum += grid->Jx_[i][j][k];
	  JySum += grid->Jy_[i][j][k];
	  JzSum += grid->Jz_[i][j][k];
	}
      }
    }

    //Calculate deposited |J|, and expected calculated |J|
    double modJ = pow(pow(JxSum,2.0)+pow(JySum,2.0)+pow(JzSum,2.0),0.5);
    double calcModJ = q/cell_volume*mod_v;

    //Run the test
    EXPECT_NEAR(modJ,calcModJ,0.00001); //Or NEAR, if calculated quantities are exact.
  }
  
}

