#include "RNG.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "assert.h"
#include <cmath>
#include <limits>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1+IMM1/RNG_NTAB)
#define RNMX (1.0-DBL_EPSILON)
#define LINE_LENGTH 2000

Random_Number_Generator::Random_Number_Generator(long int seed){
	state_ = (RNG_State*) malloc(sizeof(RNG_State));
	state_->initialSeed = seed;
  	long int* idum=&(state_->idum);
	state_->idum2=123456789;
    state_->iy=0;
    if ( *idum == 0)
		*idum=1; /* Be sure to prevent idum = 0 */
    else if(*idum < 0) 
		*idum = -(*idum);
    state_->idum2=(*idum);
    for (int j=RNG_NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      long int k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < RNG_NTAB) state_->iv[j] = *idum;
    }
    state_->iy=state_->iv[0];
	state_->generate=true;
}

Random_Number_Generator::~Random_Number_Generator(){
	free(state_);
	if(userVal_ != NULL)
		free(userVal_);
	if(userCDF_ != NULL)
		free(userCDF_);
}

double Random_Number_Generator::getUniform(){
  int j;
  long int k;
  long int* idum;
  double temp;
	
  idum=&(state_->idum);
  
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=state_->idum2/IQ2;
  state_->idum2=IA2*(state_->idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (state_->idum2 < 0) state_->idum2 += IM2;
  j=(int)(state_->iy/NDIV);              /* Will be in the range 0...RNG_NTAB-1 */
  state_->iy=state_->iv[j]-state_->idum2;                /* Here idum is shuffled, idum and idum2 */
  state_->iv[j] = *idum;                 /* are combined to generate output */
  if (state_->iy < 1) state_->iy += IMM1;
  if ((temp=AM*state_->iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NDIV
#undef RNMX

double Random_Number_Generator::getStandardNormal(){
	return getGaussian(0.0,1.0);
}

/* Adapted from Wikipedia, annotated by Denis St-Onge
 * Box Mueller generates numbers in pairs, 
 * so store both, return one at a time.
 */
double Random_Number_Generator::getGaussian(double mu, double sigma){
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	if (!state_->generate){ // just return the stored second number
		state_->generate = !state_->generate;
		return state_->z1 * sigma + mu;
	}

	double u1, u2;
	// Do loop to avoid u1 equaling to zero. 
	// Should never happen as ran2 returns between 0 and 1 exclusively
	do {
		u1 = getUniform();
		u2 = getUniform();
	} while ( u1 <= epsilon ); 

	state_->z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	state_->z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

	state_->generate = !state_->generate;

	return state_->z0 * sigma + mu;
}

void Random_Number_Generator::setRNGState(RNG_State* state){
	state_=state;
}

RNG_State* Random_Number_Generator::getRNGState(){
	return state_;
}


void Random_Number_Generator::setUserPDF(long size, double* userVal, double* userProb){
	userVal_ = (double*) calloc(size,sizeof(double));	
	userCDF_ = (double*) calloc(size,sizeof(double));	
	userSize_ = size;
	for(int i = 0; i < size; i++){
		assert(userVal[i] >= 0);
		userVal_[i] = userVal[i];
		userCDF_[i] = userProb[i];
	}
	computeUserCDF();
}

void Random_Number_Generator::loadUserPDFfromFile(char* fname){
	FILE* in = fopen(fname,"r");
	char buf[LINE_LENGTH];
	assert(in != NULL);

	userSize_=0;

	//get line count excluding comments
	fgets(buf,LINE_LENGTH,in);
	while(buf[0] != '\0' && !feof(in)){
		if(buf[0] != '#'){
			userSize_++;
		}
	}

	userVal_ = (double*) calloc(userSize_,sizeof(double));	
	userCDF_ = (double*) calloc(userSize_,sizeof(double));	

	rewind(in);

	fgets(buf,LINE_LENGTH,in);
	int i =0;
	while(buf[0] != '\0' && !feof(in)){
		if(buf[0] != '#'){
			if(i == 0)
				sscanf(buf,"%lf %lf\n",&userVal_[i],&userCDF_[i]);
			i++;
		}
	}
	computeUserCDF();
	fclose(in);
}

void Random_Number_Generator::computeUserCDF(){
	// Given PDF, calculate CDF
	userCDF_[0] = 0.0;
	for(int i = 1; i < userSize_; i++){
		userCDF_[i] += userCDF_[i-1];
	}

	// Ensure normalization (0 to 1)
	for(int i = 0; i < userSize_; i++){
		userCDF_[i] /= userCDF_[userSize_ - 1];	
	}
}

double  Random_Number_Generator::getUserNumber(){
	assert(userVal_ != NULL && userCDF_ != NULL);

	int L = 0, R = userSize_ - 1;
	int i = (R + L)/2;
	double uni = getUniform();
	//Binary Search
	while(!(userCDF_[i] <= uni && userCDF_[i + 1] > uni)){
		if(uni < userCDF_[i])
			R = i;
		else
			L = i;
		i = (R + L)/2;
	}

	//Linear Interpolation
	double len = userCDF_[i+1] - userCDF_[i];
 	double f = (uni - userCDF_[i])/len;

	return (1.0 - f)*userVal_[i] + f*userVal_[i+1];
}
