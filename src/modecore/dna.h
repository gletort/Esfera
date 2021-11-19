#ifndef __DNA_H
#define __DNA_H

#include "../tools/readXML.h"
#include "../tools/vector3d.h"
#include "../tools/random.h"

class DNA
{
	
	private:

			double condensation;
			double rate;
			Vector3d center;
			Vector3d prev_center;
			double radius;
			Random* ran;
			Vector3d cur_forces;

	public:

		DNA();
		~DNA();

		void readParameters(ReadXML* reader);

		inline void setCenter( Vector3d sph_pos ) { prev_center = center; center = sph_pos; }
		inline Vector3d* getCenter() { return &center; }
		inline void setRadius( double sph_rad ) { radius = sph_rad; }
		
		inline void resetForces() { cur_forces.set(0,0,0); }
		void addInvForce( Vector3d f );	
		inline Vector3d getForce() { return cur_forces; }

		/** \brief Test if sphere should be in contact or not */	
		double contactRate1( Vector3d spherepos, double sphererad, double time_step );
		/** \brief Test if dna from nucleole contact membrane or not */	
		double contactRateMemb( double dist, double time_step );
		/** \brief Gives Yes or No according proba of being in contact */	
		int contactProba( Vector3d spherepos );
		/** \brief Gives Yes or No according proba of being in contact if move along direction by amount howfar*/	
		int contactProbaThere( Vector3d spherepos, Vector3d dir, double howfar );
		/** \brief Gives value of proba of being in contact (between 0 and 1)*/	
		double contactProbaValue( Vector3d spherepos );
		/** \brief vector of repulsion, normalise */
		Vector3d repulsion( Vector3d spherepos);
		/** \brief If spehers is inside DNA */
		int insideDNA( Vector3d pos, double sprad);
		/** \brief attraction between a sphere and dna surface */
		Vector3d attraction( Vector3d spherepos, double rad, double cond );

		double getThickness();

		/** \brief return displacement vector of the center */
		Vector3d displacement();
};

#endif
