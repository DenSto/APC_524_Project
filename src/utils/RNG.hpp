#ifndef RNG_HPP
#define RNG_HPP

#define RNG_NTAB 32

typedef struct RNG_State{
	long int idum;
	long int idum2;
	long int iy;
	long int iv[RNG_NTAB];
	double z0,z1;
	bool generate;
}RNG_State;


class Random_Number_Generator{
	public:
		Random_Number_Generator(long seed);
		~Random_Number_Generator();
		double getUniform();
		double getStandardNormal();
		double getGaussian(double mu, double sigma);
		RNG_State *getRNGState(); 
		void setRNGState(RNG_State* state);

	private:
		RNG_State *state_;
};
#endif
