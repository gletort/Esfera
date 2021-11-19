#include "sphere_prop.h"
#include "math.h"

SphereProp::SphereProp()
{
	num_prop = 0;
	dna = NULL;
	dna_src = 0;
}

SphereProp::SphereProp(int nb)
{
	num_prop = nb;
	dna = NULL;
	dna_src = 0;
}

SphereProp::~SphereProp()
{
}

void SphereProp::readParameters(ReadXML* reader, DNA* dnaobj)
{
	// default values
	d0 = 1;
	alpha = 1;
	name = "";
	space_repul = 1;
	sph_repul = 1;
	sph_attr = 1;
	mode_diff = 0;
	max_dist = 0;
	qtotal = -1;
	qfree = 0;
	dissorate = -1;
	minvolume = 0.0042; // radius 0.1
	escape_rate = 0; // no escape
	//dna_thick = 0; // no add
	tau = -1;
	friction = 1;
	towrite = 1;

	// read parameters
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "minimal_mobility", &d0 );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "alpha", &alpha );
	(*reader).getStringValue( "sphere_prop_"+std::to_string(num_prop), "name", &name );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "membrane_repulsion", &space_repul );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "sphere_repulsion", &sph_repul );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "sphere_attraction", &sph_attr );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "maximal_distance", &max_dist );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "quantite_totale", &qtotal );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "dissolution_rate", &dissorate );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "minimal_volume", &minvolume );
	(*reader).getIntValue( "sphere_prop_"+std::to_string(num_prop), "mode_diffusion", &mode_diff );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "persistence", &tau );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "friction", &friction );
	(*reader).getIntValue( "sphere_prop_"+std::to_string(num_prop), "write", &towrite );
	
	int bpos = name.find_first_not_of(" \t");
	int epos = name.find_last_not_of(" \t");
	name = name.substr( bpos, epos );
	
	if ( qtotal > 0 )
		qfree = qtotal;	
	
	dna_stick = 0;
	dna_src = 0;
	dna_repul = sph_repul;
	dna_cond = 0;
	(*reader).getIntValue( "sphere_prop_"+std::to_string(num_prop), "dna_sticker", &dna_stick );
	(*reader).getIntValue( "sphere_prop_"+std::to_string(num_prop), "dna_source", &dna_src );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "dna_escape_rate", &escape_rate );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "dna_repulsion", &dna_repul );
	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "dna_condensation", &dna_cond );
//	(*reader).getDoubleValue( "sphere_prop_"+std::to_string(num_prop), "dna_thickness", &dna_thick );
//	if ( likedna )
	//if ( dna_src <= 0)
	dna = dnaobj;
}


void SphereProp::updateFreeQuantity( double dvol )
{
	if ( qtotal > 0 )
		qfree += dvol; 
}
