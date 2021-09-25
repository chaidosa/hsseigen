#ifndef RANDGEN_H
#define RANDGEN_H
#include<tr1/random>

#define RANDOM_GEN 
using namespace std;

typedef enum Distributions{UNIFORM_REAL, UNIFORM_INT}Distributions;

typedef tr1::mt19937                                        generator_t;
typedef tr1::uniform_real<double>                           uniform_real_distribution_t;
typedef tr1::variate_generator<generator_t, uniform_real_distribution_t> uniform_real_variate_t;

class RandGen
{
	generator_t sparfun_rand; //engine
	void*  sparfun_rand_unif; 
	void * dist; //distribution
	double GetTime();
public:
	RandGen(Distributions d);
	void SetSeed(unsigned long seed);
	void SetTimeSeed(void); 
	void SetInterval(double min, double max);
	double Rand();
	~RandGen();
};
#endif

