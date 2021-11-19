#ifndef __SPACE_H
#define __SPACE_H

#include "../tools/readXML.h"
#include "../tools/vector3d.h"

class Space
{
	
	private:

		/** \brief Main geometry radius */
		double radius;

		/** \brief other sphere radius for death star geom */
		double oradius;
		/** \brief other sphere pos for death star geom */
		Vector3d ocenter;

		/** \brief Choose geometry (0: sphere, 1: death star like */
		int geom;

		/** \brief write space contour to file */
		int write;

	public:

		Space();
		~Space();

		void readParameters(ReadXML* reader);

		/** \brief Get closest point on the space surface 
		 *
		 * Put in proj vector the vector PProj (P current position)
		 * Return if point is out of volume 2, 1 inside overlap, 0 inside no contact */
		int distToSurface( Vector3d pos, Vector3d* toSurf, double sphererad);

		/** \brief proj whe geom is a sphere */
		int distToSphere(Vector3d pos, Vector3d* toSurf, double sphererad);

		/** \brief proj whe geom is death star */
		int distToDeathStar( Vector3d pos, Vector3d* toSurf, double sphererad);

		/** \brief find closest point on anti-sphere that is inside big sphere */
		void projInsideDeathStar( Vector3d* pj, Vector3d pos );

		/** \brief distance to surface of point pt */
		double distance( Vector3d pt, bool* inside );

		/** \brief write points of the surface */
		void writeContour();
};

#endif
