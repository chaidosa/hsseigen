#include<sys/time.h>
#include"RandGen.h"

//based on https://gist.github.com/dgleich/578563

RandGen::RandGen(Distributions d)
{
	switch(d)
	{
		case UNIFORM_REAL:
				{
					sparfun_rand_unif = new uniform_real_variate_t(sparfun_rand,uniform_real_distribution_t(0.0,1.0));
				}
				break;
		default:
			printf("not supported\n");
			break;
	}
}

double RandGen::GetTime()
{
  struct timeval t; gettimeofday(&t, 0);
  return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
}

void RandGen::SetSeed(unsigned long seed)
{
  sparfun_rand.seed(seed);
  delete ((uniform_real_variate_t*)sparfun_rand_unif);
  sparfun_rand_unif = new uniform_real_variate_t(sparfun_rand, uniform_real_distribution_t(0.0, 1.0));
}

/** Return a seed based on the time. */
void RandGen::SetTimeSeed(void) 
{
    unsigned long seed = (unsigned long)GetTime();
    SetSeed(seed);
    return;
}

void RandGen::SetInterval(double min, double max)
{
  dist = new uniform_real_distribution_t(min,max); //distribution
}

double RandGen::Rand()
{
	double ret = (*((uniform_real_distribution_t *)dist))(*(uniform_real_variate_t*)sparfun_rand_unif);
	return ret;
}

RandGen::~RandGen()
{
	delete ((uniform_real_distribution_t*)dist);
	delete ((uniform_real_variate_t*)sparfun_rand_unif);
}

//TestCode
/*int main()
{
	RandGen randGen(UNIFORM_REAL);
	randGen.SetSeed(0);
	randGen.SetInterval(0.0,1.0);
	double t = randGen.Rand();
	printf("%f\n",t);
}*/
