#ifndef RNG_HPP
#define RNG_HPP

#define RNG_NTAB 32

/*!
 * Structure that contains all the infomration for a random number generator. 
 * Can be written/read using fwrite/fread.
 */
typedef struct RNG_State{
  long int initialSeed;
  long int idum;
  long int idum2;
  long int iy;
  long int iv[RNG_NTAB];
  double z0,z1;
  bool generate;
}RNG_State;


//! Class that provides methods to generate random numbers
/*!
 * 
 *  The Random Number generator class uses the ran2 algorithm from Numerical recipes.
 *  The algorithm provides fast random numbers over (0,1) exclusive with a period of over 10^15.
 * 
 *  This is then used in the implementation for numerous other distrituions (standard normal, integer...)
 *
 *
 *  This class can also be loaded with a user defined PDF, either discrete or continuous. User provided 
 *  PDF does not need to be normalized, but must be positive everywhere.
 *
 *  Discrete PDF is treat as a histogram, while the continuous PDF is treated as piecewise linear 
 *  and continous. The CDF is calculated (simple partial sum for discrete, 
 *  triangle rule for continuous) and then used for value sampling. This is done using a binary search.
 * 
 */
class Random_Number_Generator{
  public:
    Random_Number_Generator(long seed); 
    ~Random_Number_Generator();
    double getUniform();// Get a number between 0 and 1 exclusive
    long getInteger(long min, long max);// Get an integer between min and max
    double getStandardNormal(); // Get a number from a standard normal (zero mean, unit variance)
    double getGaussian(double mu, double sigma); // Get a number from an arbitrary gaussian
    RNG_State *getRNGState();  // get the current RNG state
    void setRNGState(RNG_State* state); // set the current RNG state

    void setUserPDF(bool isDiscrete,long size, double* userVal, double* userProb); // set a user-defined distribution
                                             // both userVal and userProb must be 
                                             // of size 'size'.
  
    void loadUserPDFfromFile(const bool isDiscrete, const char* fname); // load a user-defined distribution from a file
    double getUserNumber();// draw a number from user-defined distribution


  private:
    RNG_State *state_;
    double* userVal_;
    double* userCDF_;
    double* userPDF_; // needed for continous PDF for quadratic interpolation
    long userSize_;
    bool userIsDiscrete_;
    void computeUserCDF();
    void freeUserData();
};
#endif
