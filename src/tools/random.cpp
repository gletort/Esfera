#include "random.h"
#include <ctime>
#include <iostream>

Random Random::randInstance;

void Random::init()
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	/* using nano-seconds instead of seconds */
	srand((time_t)ts.tv_nsec);
	unsigned seed = rand(); //std::time(0);
	//std::cout << seed << std::endl;
	generator.seed( seed );
}

double Random::uniform()
{
	 std::uniform_real_distribution<double> dis(0.0, 1.0);
	 return dis(generator);
}

double Random::normal()
{
	 std::normal_distribution<double> dis(0.0, 1.0);
	 return dis(generator);
}

double Random::exponential()
{
	 std::exponential_distribution<double> dis(1);
	 return dis(generator);
}
