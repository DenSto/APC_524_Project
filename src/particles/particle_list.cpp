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
}

Particle_Field_List::~Particle_Field_List(){
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

    for(long ip=0;ip<np_;ip++){
        pusher_->Step(parts_[ip],parts_[ip]->field,dt);
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
