#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "../particles/particle.hpp"

class Pusher {
    public:
        virtual ~Pusher() {};
        virtual int Step(Particle* part, Field_part *field, double dt) = 0;
};



#endif
