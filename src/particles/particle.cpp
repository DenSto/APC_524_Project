#include <stdlib.h>
#include "assert.h"
#include "particle.hpp"

Particle* new_particle(){
	Particle* part = (Particle*)malloc(sizeof(Particle));
	assert(part != NULL);
	part->q=1.0;
	part->m=1.0;
	part->field=new_particle_field();
	return part;
}

void free_particle(Particle* part){
	free_particle_field(part->field);
	free(part);
	return;
}

Field_part* new_particle_field(){
	Field_part* field = (Field_part*)malloc(sizeof(Field_part));
	assert(field != NULL);
	field->e1=0.0;
	field->e2=0.0;
	field->e3=0.0;
	field->b1=0.0;
	field->b2=0.0;
	field->b3=0.0;

	return field;
}
void free_particle_field(Field_part* field){
	free(field);
	return;
}
