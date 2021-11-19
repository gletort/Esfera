#include "sphere.h"
#include "../tools/random.h"
#include <iostream>



Sphere::Sphere() : pos(0,0,0), prev_speed(0,0,0), cur_ran(0,0,0)
{
	prop = NULL;
	radius = 1;
	dna_contact = -1;
	detach = -1;
}

Sphere::Sphere( double rad, SphereProp* p, Vector3d* where, int contact ) : radius(rad), prev_speed(0,0,0), cur_ran(0,0,0), prop(p)
{
	setRadius(rad);
	pos = (*where);
	dna_contact = -1;
	detach = -1;
	prop->updateFreeQuantity(getVolume()); 

	if ( prop->dna != NULL )
	{
		dna_contact = contact;
		if ( dna_contact > 0 )
		{
			Random* rand = Random::getInstance();
			detach = rand->exponential();
			//dna_contact_point = pos;
		}
	}
}


Sphere::~Sphere()
{
}

void Sphere::loadDefault()
{
}

int Sphere::dissolve( double dt )
{
	double dvol = prop->dissorate * dt;
	if ( dvol > 0 )
	{
		// no pb
		if ( getVolume() - dvol > prop->minvolume )
		{
			removeVolume( dvol );	
			prop->updateFreeQuantity(dvol); 
			return 0;
		}
		else
		{
			// getting too small
			prop->updateFreeQuantity( getVolume() );
			//to kill
			return 1;
		}
	}
	return 0;
}

/**  Update current position from calculated velocities 
*
* Use Adam-Bashforth integration xi(t+dt) = xi(t) + 1/2*dt*(3vi(t) - vi(t-dt) ) 
* ( cf Ghaffarizadeh et al 2018, PhysiCell )
*/
void Sphere::updatePosition( double dt )
{	
	pos += 0.5 * dt	* (3*speed - prev_speed);
	prev_speed = speed;
	speed.set(0,0,0);
}


/** Current value of the parameters depending on the radius */
void Sphere::setRadius( double rad )
{ 
	radius = rad; 
	cur_diff_0 = prop->d0/ ( 6 * M_PI * radius * prop->friction );
	cur_diff_alpha = prop->alpha / ( 6 * M_PI * radius * prop->friction );
	cur_space_rep = prop->space_repul / ( 6 * M_PI * radius * prop->friction );
	//cur_sph_rep = prop->sph_repul / ( 6 * M_PI * radius );
	//cur_sph_attr = prop->sph_attr / ( 6 * M_PI * radius );
	cur_dna_repul = prop->dna_repul / ( 6 * M_PI * radius * prop->friction );
}

/** Repulsion if in contact with membrane */
void Sphere::addConfinement(Vector3d* proj, int out)
{
	// (1-|r|/R)^2
	double coef;
	if ( out == 2 )
	{		
		coef = ( 1 + proj->norm() / radius );
	}
	else
		coef = ( 1 - proj->norm() / radius );
	coef *= coef;

	speed += -1 * cur_space_rep * coef * proj->dir();
}

void Sphere::repulsion( Sphere* osphere, double dist, double deq, bool init )
{
	Vector3d dir = (osphere->pos - pos).dir();
	double coef = 1 - dist/deq;
	coef *= coef;
	//if ( init )
	//	coef *= 2;
	coef *= (prop->sph_repul+osphere->prop->sph_repul)/2.0/(6*M_PI*radius*prop->friction);
	/** in contact with DNA, repulsion less strong */
	if ( dna_contact > 0 )
	{
		if ( osphere->dna_contact > 0 )
			coef *= std::exp(-10*(prop->dna_cond));
	}
	speed -= coef * dir;
}

void Sphere::attraction( Sphere* osphere, double dist, double deq )
{
	Vector3d dir = (osphere->pos - pos).dir();
	dist -= deq;
	double coef = 1 - dist/prop->max_dist;
	coef *= coef;
	coef *= (prop->sph_attr+osphere->prop->sph_attr)/2.0/(6.0*M_PI*radius*prop->friction);
	speed +=  coef * dir;
}

void Sphere::interaction( Sphere* osphere, bool* flag_fusion, bool init )
{
	double dist = distance(osphere);
	double deq = radius + osphere->radius;

	// si proches, interaction plausible
	if ( dist < deq + prop->max_dist )
	{
		// same type
		//if ( prop == osphere->prop && dna_contact<=0 && osphere->dna_contact<=0)
		if ( prop == osphere->prop )
		{
			// cannot fuse
			if ( prop->max_dist <= 0 )
			{
				if ( dist < deq )
				{
					repulsion( osphere, dist, deq, init );
				}
			}
			else
			{
				if ( dist < deq )
				{
					if ( init )
						repulsion( osphere, dist, deq, init );
					else
						(*flag_fusion) = true;
				}
				else
				{
					if (!init )
						attraction( osphere, dist, deq );
				}
			}
		}
		// not same type, only repulsion pour le moment
		else
		{
			if ( dist < deq )
			{
				repulsion(osphere, dist, deq, init);
			}
		}
	}
}

/** Check if sphere is captured by DNA */
void Sphere::checkDNAContact(double tstep)
{
	if ( (prop->dna_stick) > 0 )
	{
	// Is attached to dna: should detach ?
	if ( detach >= 0 )
	{
		detach -= prop->escape_rate * tstep;
		// stay in contact
		if ( detach >= 0 )
		{
			// weak attraction toward dna surface
			Vector3d dnattr = (prop->dna)->attraction(pos, radius, prop->dna_cond);
			speed += 1.0/(6.0*M_PI*radius*prop->friction) * dnattr;
			//dna_contact_point += (prop->dna)->displacement();
			return;
		}
		// detach < 0, not in contact anymore
		dna_contact = 0;
	}

	// Free sphere
	/**if ( prop->dna != NULL )
	{
		// to be sure, should be always
		if ( dna_contact == 0 )
		{
			//if ( (prop->dna)->contactProba(pos) )
			if ( (prop->dna)->insideDNA(pos, radius) )
				speed += cur_dna_repul * (prop->dna)->repulsion(pos);
		}
		dna_contact = (prop->dna)->contactRate1( pos, radius, tstep );
		if ( dna_contact > 0 )
		{
			Random* rand = Random::getInstance();
			detach = rand->exponential();
		//	dna_contact_point = pos;
		}

	}*/
	}
}

void Sphere::contactNow( double distlim )
{
	if ( prop->dna_stick > 0 )
	{
		if ( dna_contact == 0 )
		{
			if ( distanceToDNA() <= distlim )
			{
				dna_contact = 1;
				detach = -1; // never detach 
			}
		}
	}
}

double Sphere::distanceToDNA()
{
	if ( dna_contact > 0 )
		return -1000;
	return (prop->dna)->distance(pos, radius);
}

/** Repuls if sphere inside DNA */
void Sphere::addRepulsionByDNA(double tstep)
{
	// not a speckle in contact
	if ( dna_contact == 0 )
	{
		double cprob = (prop->dna)->contactRate1(pos, radius, tstep); 
		if ( cprob > 0 )
		{
			Vector3d reping = cprob * cur_dna_repul * (prop->dna)->repulsion(pos);
			speed += reping;
			(prop->dna)->addInvForce(reping);
		}
	}
}

void Sphere::updateDNAForce()
{
	speed += 1.0/(6*M_PI*radius) * (prop->dna)->getForce();
}

/** Add brownian motion for given diffusion coeff */
void Sphere::addBrownianD(double dt, double activity)
{
	Vector3d ran;
	ran.randomNormal();
	// free diffusion
	if ( dna_contact <= 0 )
	{
		speed += sqrt(2*(cur_diff_0 + cur_diff_alpha * activity)/dt) * ran;
	}
	else  // diffusion of speckle in contact with dna
	{
		// small fluctuations around DNA
		speed += sqrt( 2*(cur_diff_0+cur_diff_alpha*activity)/(10*dt) ) * ran;
		// Follow DNA movement
		speed += 1.0/dt*(prop->dna)->displacement();

	}
}

/** Add DNA displacement for spheres in contact */
void Sphere::addDNADisplacement(double dt)
{
	if ( dna_contact > 0 )
	{
		// Follow DNA movement
		speed += 1.0/dt*(prop->dna)->displacement();
	}
}

/** Add persistent random motion */
void Sphere::addPersistentDiff( double dt, double activity )
{
	if ( cur_ran.norm() == 0 )
		cur_ran.randomNormal();

	Random* ran = Random::getInstance();
	if ( ran->uniform() <= dt/(prop->tau*activity) )
		cur_ran.randomNormal();

	speed += sqrt( 2*(cur_diff_0 + activity*cur_diff_alpha)/dt ) * cur_ran;
}

/** Add diffusive motion */
void Sphere::addDiffusion( double dt, double dist, double activity )
{
	switch ( prop->mode_diff )
	{
		case -1:
			// Immobile
			break;
		case 0:
			// One pure brownian
			addBrownianD( dt, activity );	
			break;
		case 1:
			// random with persistence
			addPersistentDiff( dt, activity );	
			break;
		default: 
			break;
	}	
}

/** Perform one time step 
*
* Calculate velocity from brownian motion 
* Update position
* */
void Sphere::oneStep( double dt, Space* space, double activity, bool init )
{
	Vector3d disSurf;
	int inte = space->distToSurface(pos, &disSurf, radius);
	double dis2Surf = disSurf.norm();
	
	checkDNAContact(dt);
	if ( !init )
	{
		// diffusive motion
		addDiffusion( dt, dis2Surf, activity );
	}
	else
	{
		// displace contact sphere with nucleole
		addDNADisplacement(dt);
	}

	// contact with membrane
	if ( inte )
		addConfinement(&disSurf, inte);
	else
	{
		/**if ( prop->dna_src > 0 && (dis2Surf<(3+radius)) )
		{
			double cprob = (prop->dna)->contactRateMemb( dis2Surf, dt );
			if ( cprob > 0 )
			{
				speed -= cprob * cur_dna_repul * disSurf.dir();
			}	
		}*/
	}
	
	//if ( prop->dna_src <= 0 )
	//	addRepulsionByDNA( dt );
}

void Sphere::write( std::ofstream& ofile )
{
   ofile << ";" << prop->name;
   ofile << ";" << radius;
   for ( int i = 0; i < 3; i++ )
		ofile << ";" << pos[i];	 
	ofile << ";" << dna_contact;  
}
