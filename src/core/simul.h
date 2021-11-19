/** 
 * \brief Main simulation
 * 
 * Handles initialisations
 * Reads para/init state
 */

#ifndef __SIMUL__
#define __SIMUL__

#include "../tools/random.h" 
#include "sphere_set.h" 
#include "dna.h" 
#include "../tools/readXML.h"
#include <fstream>

class Simul
{
	private:

		/** \brief simulation time step */
		double dt; 
		/** \brief simulation last time */
		double tmax; 
		/** \brief simulation initialisation time */
		double tinit; 
		/** \brief time to write results */
		double toutput; 
		/** \brief number of equilibre steps */
		int neq_steps;
		/** \brief amount of activity (fluctuations) */
		double activity;
		/** \brief time of switch (-1 no switch) */
		double tswitch;
		/** \brief amount of sphere to put in contact when switching */
		double cswitch;

		/** \brief Set of objects */
		SphereSet* all_spheres;
		/** \brief DNA object */
		DNA thedna;

		/** \brief output file */
		std::ofstream outfile;

	public:
		
		/** \brief Default constructor */
		Simul();
		/** \brief Destructor, frees memory */
		~Simul();
		
		/** \brief Run simul for all times */
		void run();

		/** \brief initialize everything */
		void initialize();

		/** \brief read parameters and/or load default values */
		void readParameters(ReadXML* reader);	
};

#endif
