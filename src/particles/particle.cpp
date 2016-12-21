#include "particle.hpp"

Particle new_particle(){
	Particle part;
	part.q=1.0;
	part.m=1.0;
	part.isGhost=0;
	part.field = new_particle_field();
	return part;
}

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
