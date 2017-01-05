#include "particle.hpp"

//! Creates a new particle struct with m=1 and everything else zeroed out 
Particle new_particle(){
	Particle part;
	part.x[0]=0.0;
	part.x[1]=0.0;
	part.x[2]=0.0;
	part.v[0]=0.0;
	part.v[1]=0.0;
	part.v[2]=0.0;
	part.gamma=0.0;
	part.q=0.0;
	part.m=1.0;
	part.isGhost=0;
	part.field = new_particle_field();
	part.my_id=0;
	part.initRank=0;
	return part;
}

//! Creates a particle-position field struct with everything zeroed out
Field_part new_particle_field(){
	Field_part field;
	field.e1=0.0;
	field.e2=0.0;
	field.e3=0.0;
	field.b1=0.0;
	field.b2=0.0;
	field.b3=0.0;

	return field;
}
