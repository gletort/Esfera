#include "simul.h"
#include <iostream>
#include <omp.h>


Simul::Simul()
{
	all_spheres = new SphereSet();
	all_spheres->setDNA( &thedna );
}

Simul::~Simul()
{
	delete all_spheres;
}

void Simul::readParameters(ReadXML* reader)
{
	(*reader).loadXMLFile( "simulation", "init_config", "space", "dna", "sphere_prop");
	// put default values
	dt = 0.1;
	tmax = 10;
	tinit = 0;
	toutput = 1;
	int nboutput = 20;
	int omp = 4;
	neq_steps = 1000;
	activity = 1.0;
	tswitch = -1.0;
	cswitch = 0.5;

	// read parameters
	(*reader).getDoubleValue( "simulation", "time_step", &dt );
	(*reader).getDoubleValue( "simulation", "final_time", &tmax );
	(*reader).getIntValue( "simulation", "nb_output", &nboutput );
	(*reader).getDoubleValue( "simulation", "initial_time", &tinit );
	(*reader).getIntValue( "simulation", "equilibre_nsteps", &neq_steps );
	(*reader).getIntValue( "simulation", "openmp", &omp );
	(*reader).getDoubleValue( "simulation", "switch_time", &tswitch );
	(*reader).getDoubleValue( "simulation", "switch_contact", &cswitch );
	
	(*reader).getDoubleValue( "simulation", "activity", &activity );

	omp_set_num_threads(omp);
	toutput = tmax/nboutput;
}

void Simul::initialize()
{
	ReadXML reader;
	Random* rand = Random::getInstance();
	rand->init();
	readParameters( &reader );
	thedna.readParameters( &reader );
	all_spheres->init( &reader, dt, neq_steps, activity );
}

void Simul::run()
{
	Random* rand = Random::getInstance();
	double t = 0;
	double tnext = 0;
	int nwrite = 0;
	while ( t < tinit )
	{
		all_spheres->step(dt, activity, false);
		t += dt;
	}
	t = 0;
	while ( t <= (tmax+dt*0.9) )
	{
		if ( tswitch > 0)
		{
			if (t>tswitch)
			{
				all_spheres->switchNow(cswitch);
				activity *= 1.8;
				tswitch=-1;
			}
		}
		// write current spheres state
		if ( (t - tnext) >= (-dt*0.99) )
		{
			tnext += toutput;
			std::string filename; 
			filename.resize( 1024 );
			sprintf( (char*) filename.c_str() , "output//spheres_%05d.txt" , (int)round(nwrite) ); 
			nwrite ++;
			std::ofstream outfile (filename.c_str(), std::ofstream::out);
			outfile << "%Time;Id;Name;Rad;X;Y;Z;DNAContact\n";
			all_spheres->write(outfile, t);
			outfile << std::endl;
			outfile.close();
		}
		
		all_spheres->step(dt, activity, false);
		t += dt;
	}

}
