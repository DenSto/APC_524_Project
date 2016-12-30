#ifndef RELA_BORIS_HPP
#define RELA_BORIS_HPP
#include "pusher.hpp"

class Relativistic_Boris : public Pusher {
    public:
        Relativistic_Boris();
        ~Relativistic_Boris();
        int Step(Particle *part, Field_part *field, double dt);
		
};

#endif
