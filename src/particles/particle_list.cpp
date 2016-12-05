#include "particle_list.hpp"
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"

Particle_Field_List::Particle_Field_List(long np){
    np_=np;

	parts_.reserve((long)1.5*np); // Have at least 1.5x the number of particles for 
	                      // slosh room. Excessive maybe?

    Field_part* fields_ = (Field_part*)malloc(sizeof(Field_part)*np_);
    assert(fields_!=NULL);
}

Particle_Field_List::~Particle_Field_List(){
    free(fields_);
}

void Particle_Field_List::Load(){
    //dummy code inserted by Yuan for testing main.cpp
    for(long ip=0;ip<np_;ip++){
		Particle* p = new_particle();
        p->x1=-1.0;
        p->x2=1.0;
        p->x3=ip*1.0;

        p->v1=0.0;
        p->v2=4.0;
        p->v3=1.0;
		parts_.push_back(p);
    } 
}

void Particle_Field_List::Push(double dt){
    Boris *boris = new Boris();

    for(long ip=0;ip<np_;ip++){
        boris->Step(parts_[ip],&(fields_[ip]),dt);
    }

}

void Particle_Field_List::Pass(){
}

long Particle_Field_List::nParticles(){
	return np_;
}

void Particle_Field_List::SortParticles(Particle_Compare comp){
	std::sort(parts_.begin(),parts_.end(),comp);
}
