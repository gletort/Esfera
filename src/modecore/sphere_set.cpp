#include "sphere_set.h"
#include <iostream>
#include "../tools/random.h"

SphereSet::SphereSet() : max_id(0)
{
	spheres.clear();
	space = NULL;
	dna_center = -1;
}

SphereSet::~SphereSet()
{
	std::map<int, Sphere*>::iterator it;
	while ( !spheres.empty() )
	{
		it = spheres.begin();
		Sphere* sph = it->second;
		spheres.erase( it->first );
		delete sph;
	}

	std::map<int, SphereProp*>::iterator itprop;
	while ( !props.empty() )
	{
		itprop = props.begin();
		SphereProp* sph = itprop->second;
		props.erase( itprop->first );
		delete sph;
	}

	if ( space )
		delete space;
}

void SphereSet::init( ReadXML* reader, double dt, int neq, double act )
{
	initSpace(reader);
	initProps( reader );
	initSpheres( reader, dt, neq, act );
}

void SphereSet::initSpace( ReadXML* reader )
{
	space = new Space();
	space->readParameters(reader);
	space->writeContour();
}

void SphereSet::addSpheresProp( ReadXML* reader, std::string pname, int np, int res, int withContact )
{
	double rmean = 0;
	reader->getDoubleValue( "init_config", "radius_mean", &rmean, pname);
	double rstd = 0;
	reader->getDoubleValue( "init_config", "radius_std", &rstd, pname);
	int mode = 0; // initial distrib of spheres
	reader->getIntValue( "init_config", "distrib_mode", &mode, pname);
	double len = 0; // of size of the distrib
	reader->getDoubleValue( "init_config", "distrib_length", &len, pname);
	double contactProp = 0; // proportion of cell in dna contact
	reader->getDoubleValue( "init_config", "dna_contact_prop", &contactProp, pname);
	int nContact = 0;
	if ( contactProp > 0 )
	{
		nContact = (int) (res * contactProp);
	}
	if ( withContact > 0 )
	{
		for ( int j = 0; j < nContact; j++ )
		{
			addOneSphere( np, rmean, rstd, mode, len, 1 );
		}
	}
	else
	{
		for ( int j = 0; j < (res-nContact); j++ )
		{
			addOneSphere( np, rmean, rstd, mode, len, 0 );
		}
	}
}


void SphereSet::addCenteredSpheresProp( ReadXML* reader, std::string pname, int np, int res )
{
	double rmean = 0;
	reader->getDoubleValue( "init_config", "radius_mean", &rmean, pname);
	double rstd = 0;
	reader->getDoubleValue( "init_config", "radius_std", &rstd, pname);
	for ( int j = 0; j < res; j++ )
	{
		addOneSphere( np, rmean, rstd, -1, 0, -1 );
	}
}


void SphereSet::replaceSpheres( ReadXML* reader, std::string pname )
{
	int mode = 0; // initial distrib of spheres
	reader->getIntValue( "init_config", "distrib_mode", &mode, pname);
	double len = 0; // of size of the distrib
	reader->getDoubleValue( "init_config", "distrib_length", &len, pname);

	Vector3d where(0,0,0);
	placeSphere( (spheres[dna_center])->getRadius(), &where, mode, len);
	(spheres[dna_center])->moveVector( where );

	if ( mode != 2) 
	{
		// replace all spheres in contact (might get push outside)
		// update
		#pragma omp parallel for 	
		for ( int i = 0; i < spheres.size(); i++ )
		{
			auto it = spheres.begin();
			std::advance(it, i);
			if ( (it->second)->inContact() )
				(it->second)->moveVector( where );
		}
	}
}

void SphereSet::initSpheres( ReadXML* reader, double dt, int neq, double activ )
{
	int res = -1;
	// do first the dna spheres, centered
	if ( dna_prop != -1 )
	{
		std::string pname = "add_"+(props[dna_prop])->getName();
		reader->getIntValue( "init_config", "number", &res, pname);
		// this type of spheres is present
		if ( res != -1 )
		{
			int mode = 0; // initial distrib of spheres
			reader->getIntValue( "init_config", "distrib_mode", &mode, pname);
			if ( mode != 2 )
				addCenteredSpheresProp( reader, pname, dna_prop, res );
			else
				addSpheresProp( reader, pname, dna_prop, res, 0 );
		}
	}

	// do other spheres, with contact
	for ( auto it : props )
	{
		res = -1;
		if ( it.first != dna_prop )
		{
			std::string pname = "add_"+((it.second)->getName());
			reader->getIntValue( "init_config", "number", &res, pname);
			// this type of spheres is present
			if ( res != -1 )
				addSpheresProp( reader, pname, it.first, res, 1 );
		}
	}

	// now place la dna sphere, et celles en contact
	if ( dna_prop != -1 )
	{
		std::string pname = "add_"+(props[dna_prop])->getName();
		int mode = 0; // initial distrib of spheres
		reader->getIntValue( "init_config", "distrib_mode", &mode, pname);
		if ( mode != 2 )
			replaceSpheres( reader, pname );	
	}

	// equilibrate positions	
	//for ( int i = 0; i < neq; i++ )
	//	step(dt/2.0, activ, true);
	
	// add spheres without contact
	for ( auto it : props )
	{
		res = -1;
		if ( it.first != dna_prop )
		{
			std::string pname = "add_"+((it.second)->getName());
			reader->getIntValue( "init_config", "number", &res, pname);
			// this type of spheres is present
			if ( res != -1 )
				addSpheresProp( reader, pname, it.first, res, 0 );
		}
	}

	for ( int i = 0; i < 2*neq; i++ )
		step(dt/2.0, activ, true);
}

void SphereSet::initProps( ReadXML* reader )
{
	int nps = reader->nbProps();
	dna_prop = -1;

	for ( int nprop = 0; nprop < nps; nprop ++ )
	{
		SphereProp* sp = new SphereProp(nprop);
		sp->readParameters(reader, mydna);
		props[nprop] = sp;
		if ( (props[nprop])->isDNASource() )
			dna_prop = nprop;
	}
}

void SphereSet::addOneSphere( int nprop, double rmean, double rstd, int mode, double len, int contact )
{
	Random* ran = Random::getInstance();
	double rad = 0;
	// to be sure it doesn't get <0
	while ( rad <= 0.1 )
		rad = rmean + sqrt(rstd) * ran->normal();
	Vector3d where(0,0,0);
	int ntry = 0;
	int inside = 0;
	int cont = -1;
	int maxtry = 5000;
	if ( contact != -1 )
		cont = !contact;
	while (  !inside || (cont != contact) )
	{	
		placeSphere( rad, &where, mode, len);
		inside = isInsideSpace( rad, where );
		if ( contact != -1 )
		{
			if ( mode != 2 )
				cont = mydna->contactProba( where );
			else
				cont = contact;
		}
		
		if ( inside && (cont==contact) )
		{
			// if touch other sphere, but already tried a lot, let it like that
			// otherwise try to replace 
			if ( touchOtherSphere( where, rad ) )
			{
				if ( ntry <= maxtry*0.95 )
					inside = 0;
			}
		}

		ntry++;
		if ( ntry >= maxtry )
			break;
	}
	spheres[max_id] = new Sphere(rad, (props[nprop]), &where, cont );
	max_id ++;

	// update DNA infos
	if ( dna_center < 0 )
	{
		if ( (props[nprop])->isDNASource() )
		{
			dna_center = max_id-1;
			mydna->setCenter( where );
			mydna->setRadius( rad );
		}
	}
}

int SphereSet::isInsideSpace( double rad, Vector3d where )
{
	// is inside Space
	Vector3d tmp(0,0,0);
	int res = space->distToSurface( where, &tmp, rad );
	if ( res >= 2 )
		return 0;
	return 1;
}	

void SphereSet::placeSphere( double rad, Vector3d* where, int mode, double length )
{
	Random* ran = Random::getInstance();
	double rlim = length-rad;
	switch(mode)
	{
		// initial centered
		case -1:
			where->randomAndNormed();
			(*where) = 0*(*where);
			break;
		// sphere center in O of radius length
		case 0:
			where->randomAndNormed();
			(*where) = pow(ran->uniform(),1.0/3.0)*rlim*(*where);
			break;
		// distance equi-proba (not position)
		case 1:
			where->randomAndNormed();
			(*where) = ran->uniform()*rlim*(*where);
			break;
		// sphere center in O of radius length, independant of contact
		case 2:
			where->randomAndNormed();
			(*where) = pow(ran->uniform(),1.0/3.0)*rlim*(*where);
		default:
			break;
	}
}

/** Check if overlap */
bool SphereSet::touchOtherSphere( Vector3d where, double rad )
{
	for ( auto sph : spheres )
	{
		if ( sph.second->distanceVec(where) < (rad+sph.second->getRadius())  )
			return true;
	}
	return false;
}

/** perform one time step for each sphere */
void SphereSet::step( double dt, double activity, bool initial )
{
	bool fusion;
	std::map<int, int> tofuse;
	tofuse.clear();
	
	// update DNA center infos
	if ( dna_center >= 0 )
	{	
		// update
		mydna->setCenter( spheres[dna_center]->getPosition() );
		mydna->setRadius( spheres[dna_center]->getRadius() );
		mydna->resetForces();
	}

	// individual update
	#pragma omp parallel for	
	for ( int i = 0; i < spheres.size(); i++ )
	{
		auto it = spheres.begin();
		std::advance(it, i);
		// find if contact with other spheres
		for ( auto osph : spheres )
		{
			// remove self case
			if ( it->first != osph.first )
			{
				bool fusion = false;
				(it->second)->interaction(osph.second, &fusion, initial);
				if ( !initial && fusion )
				{
					#pragma omp critical (fusing)
					{// not already taken
					if ( (tofuse.count( it->first ) <= 0) && (tofuse.count(osph.first) <=0) )
					{
						tofuse[it->first] = osph.first;
						tofuse[osph.first] = it->first;
					}
					}
				}
			}
		}
		(it->second)->oneStep(dt, space, activity, initial);
	}

	// add calculated spheres-dna interaction		
	// spheres[dna_center]->updateDNAForce();

	// update
	#pragma omp parallel for	
	for ( int i = 0; i < spheres.size(); i++ )
	{
		auto it = spheres.begin();
		std::advance(it, i);
		(it->second)->updatePosition(dt);
	}

	if ( !initial ) 
		doFusion( tofuse );
	//dissolveStep(dt);
}

void SphereSet::doFusion( std::map<int, int> fusing )
{
	// Do fusions now if needed
	while ( !fusing.empty() )
	{
		int sphere1 = fusing.begin()->first;
		int sphere2 = fusing.begin()->second;
		Sphere* sph1 = spheres[sphere1];
		Sphere* sph2 = spheres[sphere2];
		double rad1 = sph1->getRadius();
		double rad2 = sph2->getRadius();

		// merge sp2 into sph1
		if ( rad1 > rad2 )
		{
			sph1->setRadius( pow(rad1*rad1*rad1+rad2*rad2*rad2, 1.0/3.0) );
			spheres.erase(sphere2);
			delete sph2;			
		}
		// merge sp1 into sp2
		else
		{
			sph2->setRadius( pow(rad1*rad1*rad1+rad2*rad2*rad2, 1.0/3.0) );
			spheres.erase(sphere1);
			delete sph1;			
		}

		fusing.erase(sphere1);
		fusing.erase(sphere2);
	}
}

void SphereSet::dissolveStep(double dt)
{
	std::vector<int> tokill;
	tokill.clear();
	int killing;

	// dissolve each sphere according to its own rate
	for ( auto sph : spheres )
	{
		killing = (sph.second)->dissolve(dt);
		if ( killing )
			tokill.push_back( (sph.first) );
	}

	// remove totally dissolved cells
	for ( int i = 0; i < tokill.size(); i ++ )
	{
		Sphere* sph = spheres[ tokill[i] ];
		spheres.erase(tokill[i]);
		delete sph;
	}
}

void SphereSet::write( std::ofstream& ofile, double tcur )
{
	bool first = true;
	for ( auto sph: spheres )
	{
		if ( (sph.second)->toWrite() > 0 )
		{
			if (!first) ofile << "\n";
			else first = false;
			ofile << tcur;
			ofile << ";" << (sph.first);
			(sph.second)->write( ofile );
		}
	}
}
