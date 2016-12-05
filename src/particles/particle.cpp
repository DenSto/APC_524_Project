#include <stdlib.h>
#include<assert.h>
#include "particle.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"

Particle* new_particle(){
	Particle* part = (Particle*)malloc(sizeof(Particle));
	part->q=1.0;
	part->m=1.0;
	return part;
}

void free_particle(Particle* part){
	free(part);
	return;
}

Field_part* new_particle_field(){
	Field_part* field = (Field_part*)malloc(sizeof(Field_part));
	return field;
}
void free_particle_field(Field_part* field){
	free(field);
	return;
}

Particle_Field_List::Particle_Field_List(int np){
    np_=np;

    Particle *parts_ = (Particle*)malloc(sizeof(Particle)*np_);
    assert(parts_!=NULL);

    Field_part* fields_ = (Field_part*)malloc(sizeof(Field_part)*np_);
    assert(fields_!=NULL);
}

Particle_Field_List::~Particle_Field_List(){
    free(parts_);
    free(fields_);
}

void Particle_Field_List::Load(){
    //dummy code inserted by Yuan for testing main.cpp
    int ip;
    for(ip=0;ip<np_;ip++){
        parts_[ip].x1=-1.0;
        parts_[ip].x2=1.0;
        parts_[ip].x3=ip*1.0;

        parts_[ip].v1=0.0;
        parts_[ip].v1=0.0;
        parts_[ip].v3=1.0;
    } 
}

void Particle_Field_List::Push(double dt){
    Boris boris;
    int ip;
    for(ip=0;ip<np_;ip++){
        boris.Step(&(parts_[ip]),&(fields_[ip]),dt);
    }

}

void Particle_Field_List::Pass(){
}
