#include "dna.h"

DNA::DNA()
{
	ran = Random::getInstance();
}

DNA::~DNA()
{
}

void DNA::readParameters( ReadXML* reader )
{

	//default values
	condensation = 3;
	rate = 1;
	// read
	reader->getDoubleValue("dna", "condensation", &condensation);
	reader->getDoubleValue("dna", "rate", &rate);

}

void DNA::addInvForce( Vector3d f )
{
	#pragma omp critical (dnarep )
	{
		cur_forces -= f;
	}
}

/** Test if sphere should be in contact or not */	
double DNA::contactRate1( Vector3d spherepos, double sphererad, double time_step )
{
	// not too be sensitive on time step
	//if ( ran->uniform() <= time_step/10 ) 
	{
		double dist = (spherepos - center).norm() - sphererad;
		if ( dist <= radius ) return 1;  // inside nucleole
		double func = rate * time_step * std::exp(-(dist-radius)/condensation);
		//double func = std::exp(-(dist-radius)/condensation);
		//std::cout << "try " << dist << " " << func << std::endl;
		if ( ran->uniform() <= func )
			return func;
	}
	return 0;
}

/** Test if dna from nucleole contact membrane or not */	
double DNA::contactRateMemb( double dist, double time_step )
{
	// not too be sensitive on time step
	//if ( ran->uniform() <= time_step/10 ) 
	{
		double func = std::exp(-(dist-radius)/condensation);
		//std::cout << "try " << dist << " " << func << std::endl;
		if ( ran->uniform() <= func )
			return func;
	}
	return 0;
}

/** Test if sphere should be in contact or not */	
int DNA::contactProba( Vector3d spherepos )
{
	double dist = (spherepos - center).norm();
	double func = std::exp(-(dist-radius)/condensation);
		//std::cout << "try " << dist << " " << func << std::endl;
	return ( ran->uniform() <= func );	
}

double DNA::getThickness()
{
	double dist = ran->uniform()*radius*0.5;
	int ntry = 0;
	while ( ran->uniform() > (std::exp(-(dist-radius)/condensation)) )
	{
		dist = ran->uniform()*radius*0.5;
		ntry++;
		if (ntry > 100) return 0;
	}
	return dist;
}

/** Test if sphere will be in contact or not if go in vec direction */	
int DNA::contactProbaThere( Vector3d spherepos, Vector3d dir, double howfar )
{
	spherepos = spherepos + howfar*dir;
	double dist = (spherepos - center).norm();
	if ( dist >= radius )
	{
		double func = std::exp(-(dist-radius)/condensation);
		//std::cout << "try " << dist << " " << func << std::endl;
		return ( ran->uniform() <= func );	
	}
	else 
		return 1;
}

/** Return proba of contact */
double DNA::contactProbaValue( Vector3d spherepos )
{
	double dist = (spherepos - center).norm();
	if ( dist >= radius )
	{
		double func = std::exp(-(dist-radius)/condensation);
		//std::cout << "try " << dist << " " << func << std::endl;
		return ( func );	
	}
	return 1;
}

int DNA::insideDNA( Vector3d pos, double sprad )
{
	double dist = (pos - center).norm();
	return ( (dist-sprad) <= radius);
}

Vector3d DNA::repulsion( Vector3d spherepos )
{
	Vector3d diff = spherepos - center;
	//double dist = diff.norm();
	diff.normalize();
	//diff *= std::exp(-(dist-radius)/condensation);
	return diff;
}

Vector3d DNA::attraction( Vector3d spherepos, double rad, double cond )
{
	double deq = rad + radius;
	Vector3d diff = (spherepos - center);
	// if not too close
	if ( diff.norm() >= deq )
	{
		double dist = diff.norm() - deq;
		double coef = 1 - 1.0/(1.0+dist);  //0 if 0 and tend to 1 if infty
		return cond*coef*(diff.dir());
	}
	return (0*diff.dir()); // small additionnal repulsion
}

/** return displacement vector of the center */
Vector3d DNA::displacement()
{
	return (center-prev_center);
}
