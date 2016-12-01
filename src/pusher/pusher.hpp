#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "../particles/particle.hpp"

class Pusher {
    public:
        Pusher();
        ~Pusher();
        virtual int Step(Particle* part, Field_part *field, double dt);
};



#endif
