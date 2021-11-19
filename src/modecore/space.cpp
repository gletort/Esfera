#include "space.h"
#include <fstream>
#include <sstream>

Space::Space() : ocenter(0,0,0)
{
}

Space::~Space()
{
}

void Space::readParameters( ReadXML* reader )
{

	//default values
	geom = 0;
	radius = 10;
	oradius = 5;
	double oycenter = -12;
	write = 0;

	// read
	reader->getIntValue("space", "geometry", &geom);
	reader->getDoubleValue("space", "radius", &radius);
	if ( geom == 1 )
	{
		reader->getDoubleValue("space", "deathstar_radius", &oradius);
		reader->getDoubleValue("space", "deathstar_ypos", &oycenter);
	}
	reader->getIntValue("space", "write", &write);

	ocenter[0] = 0;
	ocenter[1] = oycenter;
	ocenter[2] = 0;
}
		
/** Get closest point on the space surface 
*
* Put in proj vector the vector PProj (P current position)
* return distance to surface mode: 0, inside no contact, 1 inside overlap, 2 outside */
int Space::distToSurface( Vector3d pos, Vector3d* toSurf, double sphererad  )
{
	switch (geom)
	{
		case 0 :
			return distToSphere( pos, toSurf, sphererad );
			break;
		case 1 :
			return distToDeathStar( pos, toSurf, sphererad );
			break;
		default:
			return distToSphere( pos, toSurf, sphererad);
			break;
	}
}

int Space::distToSphere( Vector3d pos, Vector3d* toSurf, double sphererad )
{
	double norm = pos.norm();
	(*toSurf) = pos.dir();
	(*toSurf) *= (radius-norm);
	// is even out of the volume
	if ( norm > radius )
	{
		(*toSurf) *= -1;
		return 2; // outside
	}
	return ( (norm+sphererad) >= radius ); // 1 if contact, 0 otherwise
}


int Space::distToDeathStar( Vector3d pos, Vector3d* toSurf, double sphererad )
{
	double outsidemain = pos.norm() - radius;
	double outsideanti = (pos-ocenter).norm() - oradius;

	// inside anti, proj on anti surface
	if ( outsideanti < 0 )
	{
		Vector3d proj;
		projInsideDeathStar(&proj, pos);
		(*toSurf) = 1*(pos-proj);
		return 2;
	}

	// outside anti, inside main	
	if ( outsidemain < 0 )
	{
		// closest win
		// plus proche et inside de la main sphere
		if ( (-outsidemain) < outsideanti  )
		{
			(*toSurf) = -outsidemain * pos.dir();
			return (pos.norm()+sphererad) >= radius;
		}
		
		// plus proche et outside anti
		(*toSurf) = -outsideanti * (pos-ocenter).dir();
		return ( outsideanti <= sphererad ); 
	}

	// outside anti, outside main
	(*toSurf) = outsidemain * pos.dir();
	return 2;
	
}

/** find closest point on anti-sphere that is inside big sphere */
void Space::projInsideDeathStar( Vector3d* pj, Vector3d pos )
{
	int Ndiscr = 100;
	double step = 2.0 * M_PI / Ndiscr;
	double theta, phi;
	double x,y,z;
	double mdist = 10000000;
	Vector3d proj(0,0,0);
	// parcourt toute la surface de l'anti-sphere
	for ( int itheta = 0; itheta < Ndiscr; itheta++ )
	{
		for ( int iphi = 0; iphi < Ndiscr/2; iphi ++ )
		{
				theta = itheta * step;
				phi = iphi * step;
				x = oradius * sin(phi) * cos(theta);
				y = ocenter[1] + oradius * sin(phi) * sin(theta);
				z = oradius * cos(phi);
				
				proj.set(x, y, z);
				// point inside, define contour
				if ( proj.norm() <= radius )
				{
					double dist = (proj-pos).norm();
					if ( dist <= mdist )
					{
						mdist = dist;
						(*pj) = proj; 
					}
				}
			}
	}
}


double Space::distance( Vector3d pt, bool* inside )
{
	Vector3d proj;
	(*inside) = (distToSurface( pt, &proj, 0) < 2);
	return proj.norm();
}

/** write points of the surface */
void Space::writeContour()
{
	if ( !write )
		return;

	// write space contour
	std::string filename("output/space.txt"); 
	std::ofstream outfile (filename.c_str(), std::ofstream::out);
	std::ios::sync_with_stdio(false);
	outfile << "%X;Y;Z\n";
	std::stringstream cur;
	double lim = radius * 1.15; // in case more complicated geom later
	int Ndiscr = 100;
	double step = 2.0* lim / Ndiscr;
	double x,y,z;
	double inx = 0;
	double iny = 0;
	double inz = 0;
	int npt = 0;
	double d;
	bool in;
	for ( int ix = 0; ix < Ndiscr; ix++ )
	{
		for ( int iy = 0; iy < Ndiscr; iy ++ )
		{
			cur.str(std::string());
			for ( int iz = 0; iz < Ndiscr; iz++ )
			{
				x = -lim + ix*step;
				y = -lim + iy*step;
				z = -lim + iz*step;
				
				d = distance( Vector3d(x,y,z), &in );
				if ( in )
				{
					inx += x;
					iny += y;
					inz += z;
					npt ++;
				}
				if ( d <= 0.025 )
					cur << x << ";" << y << ";" << z << "\n";
			}
			outfile << cur.str();
		}
	}	
	
	outfile.close();
	std::cout << "********************************************** " << std::endl;
	std::cout << "Space file written " << std::endl;
	std::cout << "Centroid pos " << "\t" << (inx/npt) << "\t" << (iny/npt) << "\t" << (inz/npt) << std::endl;
	std::cout << "Vol " << "\t" << (npt*step*step*step) << std::endl;
	std::cout << "********************************************** " << std::endl;
}
