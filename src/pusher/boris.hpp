#ifndef BORIS_HPP
#define BORIS_HPP
#include "pusher.hpp"

class Boris : public Pusher {
    public:
        Boris();
        ~Boris();
        int Step(Particle *part, Field_part *field, double dt);
		
};

#endif
