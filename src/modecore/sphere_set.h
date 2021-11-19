#ifndef __SPHERE_SET__
#define __SPHERE_SET__

#include <map>
#include "sphere.h"
#include "space.h"
#include "../tools/readXML.h"

class SphereSet
{
	private:
		/** \brief map of objects, classed by id number */
		std::map<int, Sphere*> spheres;
		/** \brief vector of  */
		std::map<int, SphereProp*> props;
		/** \brief max id used so far */
		int max_id;
		/** \brief Space of the simulation */
		Space* space;
		/** \brief DNA of the simulation */
		DNA* mydna;		
		/** \brief number of the dna center sphere (nucleole) */
		int dna_center;
		/** \brief number of the dna sphere property */
		int dna_prop;

	public:

		SphereSet();
		~SphereSet();

		/** \brief set dna pointer */
		inline void setDNA( DNA* d ) { mydna = d; }
		/** \brief initialize everything */
		void init( ReadXML* reader, double dt, int neq, double acti );
		/** \brief initialize space */
		void initSpace(ReadXML* reader);
		/** \brief initialize all the properties */
		void initProps(ReadXML* reader);
		/** \brief initialize spheres from the para file */
		void initSpheres(ReadXML* reader, double dt, int neq, double activ);
		/** \brief Add all spheres of one given type */
		void addSpheresProp( ReadXML* reader, std::string pname, int np, int res, int withContact );
		void addCenteredSpheresProp( ReadXML* reader, std::string pname, int np, int res );
		void replaceSpheres( ReadXML* reader, std::string pname );
		/** \brief add one sphere and its id 
		 *
		 * @nprop number of the property of new sphere
		 * @rmin, rmax: get random radius in [rmin, rmax], uniform distrib */
		void addOneSphere( int nprop, double rmean, double rstd, int mode, double length, int contact );
		/** \brief Find where to place the new sphere
		 *
		 * Mode: 0 inside central sphere of radius length */
		void placeSphere( double rad, Vector3d* where, int mode, double length );
		
		/** \brief if Position where with radius rad is inside space (or in contact) */
		int isInsideSpace( double rad, Vector3d where );
		
		/** \brief Return if given position overlap with already placed spheres */
		bool touchOtherSphere( Vector3d where, double rad );

		/** \brief perform one time step for each sphere */
		void step( double dt, double activity, bool init );
		/** \brief fuse spheres if needed */
		void doFusion( std::map<int, int> fusing );
		/** \brief Dissolve individual spheres if needed, remove too small spheres */
		void dissolveStep( double dt );

		/** \brief add results to file stream */
		void write( std::ofstream& ofile, double tcur );

};

#endif
