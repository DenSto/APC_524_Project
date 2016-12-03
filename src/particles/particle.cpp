#include <stdlib.h>
#include "particle.hpp"


Particle* new_particle(){
	Particle* part = malloc(sizeof(Particle));
	part->q=1.0;
	part->m=1.0;
	return part;
}

void free_particle(Particle* part){
	free(part);
	return;
}

Field_part* new_particle_field(){
	Field_part* ret = malloc(sizeof(Field_part));
	return ret;
}
void free_particle_field(Field_part* field){
	free(field);
	return;
}
