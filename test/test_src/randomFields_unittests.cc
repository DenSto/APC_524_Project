#define MAIN_CPP
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gtest/gtest.h"
#include "math.h"
#include "RNG.hpp"
#include "../../src/grid/grid.hpp"

#include "field_bc_factory.hpp"

#include "particle_handler.hpp"

#define SQR(x) ((x)*(x))

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// tests requiring internal examination
class RandomFieldTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    //MPI
    rank_MPI = 0;
    size_MPI = 1;
    //Read input file
	
	dt =0.01;

    Input *input =  new Input();
    char filename[200];
    sprintf(filename, "test_data/randomField_unittest_data.txt");
    input->readinfo(filename);
    input_info = input->getinfo();

    //Initialize the domain, grid, and particle_handler
    domain = new Domain(input_info);
    grid = new Grid(domain->getnxyz(),1,domain->getxyz0(),domain->getLxyz());
    part_handler = new Particle_Handler();
    part_handler->Load(input_info,domain);

	grid->constJ(0,0,0);
	grid->constE(0,0,0);
	grid->constB(0,0,0);
	grid->constRho(0);

	srand(time(NULL));

	RNG = new Random_Number_Generator(rand());
  }

  virtual void TearDown() {
    delete grid;
	delete domain;
	delete part_handler;
	delete RNG;
  }

  double dt;
  Grid *grid;
  Random_Number_Generator* RNG;
  Particle_Handler *part_handler;
  Domain *domain;
  Input_Info_t *input_info;
};

// Test field interpolation is working.
TEST_F(RandomFieldTest, testEnergy) {

	double**** fields = grid->getFieldPtr();
	int ind[6],beg[6],size[3];
	for(int n = 0; n < 6; n++){
		grid->getRealIndices(n,ind);
		for(int i = ind[0]; i < ind[3]; i++){
			for(int j = ind[1]; j < ind[4]; j++){
				for(int k = ind[2]; k < ind[5]; k++){
					fields[n][i][j][k] = RNG->getStandardNormal();
				}	
			}
		}
	} 

	Field_BC_Factory::getInstance().Construct(domain,grid,input_info);

	int nt = input_info->nt;
	FILE* fp = fopen("testit","w");
	grid->getRealIndices(0,beg);
	grid->getnxyzPhys(size);
	for(int n = 0; n < nt; n++){
		double EEx = 0;
		double EEy = 0;
		double EEz = 0;
		double EBx = 0;
		double EBy = 0;
		double EBz = 0;
		double vol = grid->getCellVolume();

		for(int i = beg[0]; i < beg[0] + size[0]; i++){
			for(int j = beg[1]; j < beg[1] + size[1]; j++){
				for(int k = beg[2]; k < beg[2] + size[2]; k++){
					EEx += 0.5*vol*SQR(0.5*(fields[E_X][i][j][k] + fields[E_X][i+1][j][k]));
					EEy += 0.5*vol*SQR(0.5*(fields[E_Y][i][j][k] + fields[E_Y][i][j+1][k]));
					EEz += 0.5*vol*SQR(0.5*(fields[E_Z][i][j][k] + fields[E_Z][i][j][k+1]));
					EBx += 0.5*vol*SQR(0.5*(fields[B_X][i][j][k] + fields[B_X][i+1][j][k]));
					EBy += 0.5*vol*SQR(0.5*(fields[B_Y][i][j][k] + fields[B_Y][i][j+1][k]));
					EBz += 0.5*vol*SQR(0.5*(fields[B_Z][i][j][k] + fields[B_Z][i][j][k+1]));
				}
			}
		}
		
		double B_en = EBx + EBy + EBz;
		double E_en = EEx + EEy + EEz;
		double E_tot = B_en + E_en;

		double EE[6];

		for(int q = 0; q < 6; q++){
			EE[q]=0;
			grid->getRealIndices(q,ind);
			for(int i = ind[0]; i < ind[3]; i++){
				for(int j = ind[1]; j < ind[4]; j++){
					for(int k = ind[2]; k < ind[5]; k++){
						EE[q] += 0.5*vol*SQR(fields[q][i][j][k]);
					}
				}
			}
		}
		double B_en1 = EE[3] + EE[4] + EE[5];
		double E_en1 = EE[0] + EE[1] + EE[2];
		double E_tot1 = B_en1 + E_en1;
		fprintf(fp,"%d %f %f %f %f %f %f\n",n,B_en,E_en,E_tot, B_en1, E_en1, E_tot1);

		grid->executeBC(-2);	

		grid->evolveFields(dt);
		grid->executeBC(-1);
	}
	fclose(fp);
}

#undef SQR
