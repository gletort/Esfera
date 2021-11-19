/**
 * \brief diffusive spheres
 *
 * Has a radius, can diffuse
 */

#ifndef __SPHERE__
#define __SPHERE__

#include "../tools/vector3d.h"
#include "sphere_prop.h"
#include <fstream>
#include "space.h"

class Sphere
{
	protected:
		/** \brief sphere current radius */
		double radius;
		/** \brief current position */
		Vector3d pos;
		/** \brief previous velocity */
		Vector3d prev_speed;
		/** \brief current velocity */
		Vector3d speed;
		/** \brief indicates if is in contact/bloqued with dna or not */
		int dna_contact;
		/** \brief time left to detach from dna when captured */
		double detach;
		/** \brief when in contact with DNA, position of point of contact */
		//Vector3d dna_contact_point;	
		/** \brief Indicates if can fuse with other spheres or not in current state*/
		//int lets_fuse;
		/** \brief Common parameters */
		SphereProp* prop;
		/** \brief Current value of sphere diffusion indep of activ */
		double cur_diff_0;
		/** \brief Current value of sphere diffusion factor of activ */
		double cur_diff_alpha;
		/** \brief Current value of sphere space repulsion parameter */
		double cur_space_rep;
		/** \brief Current value of sphere-sphere repulsion parameter */
		//double cur_sph_rep;
		/** \brief Current value of sphere-sphere attraction parameter */
		//double cur_sph_attr;
		/** \brief Current value of sphere-dna repul parameter */
		double cur_dna_repul;
		/** \brief Current random vector, direction of persistent Brown */
		Vector3d cur_ran;

	public:

		Sphere();
		Sphere(double rad, SphereProp* p, Vector3d* where, int contact);
	    ~Sphere();
	
		/** \brief return current sphere position */	
		inline Vector3d getPosition() { return pos; };
		/** \brief Set value of radius and corresponding parameters */
		void setRadius(double rad);
		/** \brief return radius value */
		inline double getRadius() { return radius; };
		/** \brief return current volume */
		inline double getVolume() { return (4.0/3.0* M_PI * radius *radius * radius); };
		/** \brief Decrease radius by amount for vol */
		inline double removeVolume(double vol) 
		{ setRadius(radius-pow(vol*3.0/4.0/M_PI, 1.0/3.0)); };	
		inline bool inContact()
		{ return dna_contact==1; };

		inline int toWrite()
		{ return prop->towrite; };

		inline void moveVector( Vector3d vec ) 
		{ pos += vec; };

		void updateDNAForce();

		/** \brief dissolve sphere radius by given amount by time step for dyna spheres */
		int dissolve(double dt); 

		/** \brief read parameters and initialize values */
		void loadDefault();

		/** \brief Update current position from calculated velocities 
		 *
		 * Use Adam-Bashforth integration xi(t+dt) = xi(t) + 1/2*dt*(3vi(t) - vi(t-dt) ) 
		 * ( cf Ghaffarizadeh et al 2018, PhysiCell )
		 */
		void updatePosition( double dt );

		inline double distance(Sphere* sph) { return (sph->pos - pos).norm() ; }
		inline double distanceVec(Vector3d vec) { return (vec - pos).norm() ; }

		/** \brief add repulsion */
		void repulsion( Sphere* osphere, double dist, double deq, bool init = false );
		/** \brief add attraction */
		void attraction( Sphere* osphere, double dist, double deq );

		/** \brief Add interaction (repulsion) with other sphere if close enough */
		void interaction( Sphere* osphere, bool* flag_fusion, bool init );
		
		/** \brief Check if sphere is captured by DNA */
		void checkDNAContact(double tstep);

		/** \brief Put sphere in contact with proba contact if dna sticker and not attached */
		void contactNow( double distlim );

		/** \brief Distance of sphere surface to DNA */
		double distanceToDNA();

		/** \brief Repuls if sphere inside DNA */
		void addRepulsionByDNA(double tstep);

		/** \brief Add DNA displacement for spheres in contact */
		void addDNADisplacement(double dt);

		/** \brief Repulsion if in contact with membrane */
		void addConfinement(Vector3d* disSurf, int out);

		/** \brief Perform one time step 
		 *
		 * Calculate velocity from brownian motion 
		 * */
		void oneStep( double dt, Space* space, double activity, bool init );

		/** \brief Add diffusive part of motion */
		void addDiffusion( double dt, double dist, double activity );
		
		/** \brief Add brownian motion for given diffusion coeff */
		void addBrownianD(double dt, double activity);
		
		/** \brief Add persistent motion for given diffusion coeff */
		void addPersistentDiff(double dt, double activity);

		/** \brief write spheres characteristics to file */
		void write( std::ofstream& ofile );
};

#endif
