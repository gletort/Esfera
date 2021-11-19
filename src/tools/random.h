#ifndef __RANDOM__
#define __RANDOM__

#include <random>

class Random
{
	private:

		std::mt19937 generator;
		static Random randInstance;

	public:

		Random() {};
		~Random() {};

		/** \brief seed the generator */
		void init();

		/** \brief return instance */
		static Random* getInstance() { return &randInstance; }

		/** \brief Return a double between 0 and 1 with uniform distribution */
		double uniform();
		/** \brief Return a double with normal distribution N(0,1) */
		double normal();
		/** \brief Return a double with exponential distribution exp(-x) */
		double exponential();
};

#endif
