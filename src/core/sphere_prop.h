/**
 * \brief common properties of spheres
 *
 * Has a radius, can diffuse
 */

#ifndef __SPHERE_PROP__
#define __SPHERE_PROP__

#include <string>
#include "../tools/readXML.h"
#include "dna.h"

class SphereProp
{
	friend class Sphere;

	protected:

		/** \brief Number of the property in the prop lists */
		int num_prop;	
		/** \brief Name of the sphere type (e.g. Nucleole) */
		std::string name;

		/** \brief minimal diff coefficient: KbT/eta */
		double d0;
		/** \brief coeff of damping */
		double alpha;
		/** \brief repulsion coefficient contact with membrane */
		double space_repul;
		/** \brief repulsion coefficient contact with other spheres */
		double sph_repul;
		/** \brief attraction coefficient with other spheres */
		double sph_attr;
		/** \brief Mode of diffusion: 0 uniform, */
		int mode_diff;
		/** \brief Maximal distance for spheres to interact */
		double max_dist;
		/** \brief resistance to mvt higher than standard */
		double friction;
		/** \brief quantity totale of speckles subunits 
		 *
		 * Default -1, and if stays -1 doesn't use the counter*/
		double qtotal;
		/** \brief current free quantity of speckles subunits */
		double qfree;
		/** \brief Dissolution rate if relevant */
		double dissorate;
		/** \brief Minimum volume to keep a sphere */
		double minvolume;
		
		/** \brief is nucleole (dna around it) */
		int dna_src;
		/** \brief rate of escaping from dna link */
		double escape_rate;
		/** \brief can bind to DNA */
		int dna_stick;
		/** \brief if feels contact with dna, pointer to it */
		DNA* dna;
		/** \brief repulsion coefficient contact with dna */
		double dna_repul;
		/** \brief For nucleole, add thickness due to DNA presence */
		//double dna_thick;
		/** \brief If in contact with DNA, condensate coefficient, 0 don't condense, 1 strongly */
		double dna_cond;

		/** \brief Persistence time */
		double tau;
		/** \brief If write position in output files */
		int towrite;

	public:

		SphereProp();
		SphereProp(int nb);
	    ~SphereProp();

		inline std::string getName() { return name; }

		void readParameters(ReadXML* reader, DNA* dnaobj);
		//void writeParameters();
		
		inline int isDNASource() { return dna_src; }
		
		/** \brief Update the quantity of free subunits, in Volume */
		void updateFreeQuantity( double dvol );
};

#endif
