#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>

typedef struct Particle {
	double x1,x2,x3;  /* coordinate in X,Y,Z */	
	double v1,v2,v3;  /* velocity in X,Y,Z */	
	
	double q;		  /* charge */
	double m;		  /* mass */
	
	int my_id; 		  /* particle id */
	short isGhost;    /* is particle in ghost cell? */

} Particle;

typedef struct Field_part {
    double e1, e2, e3; /* Electric field components */
    double b1, b2, b3; /* Magnetic field components */
} Field_part;

class Particle_Field_List {
    public:
        Particle_Field_List(int np); // list of np particles and their fields
        ~Particle_Field_List();
        void Load();      // Initialize particles
        void Push(double dt); // Push all particles
        void Pass();          // Pass particles accross MPI boundary
		void nParticles();

		void sortParticles(); // quicksort particle list

    private:
        int np_;
		std::vector<Particle*> parts_;    /* Vector of particles */
        Field_part *fields_; /* Field values at particle locations */
};

Particle* new_particle();
void free_particle(Particle* part);

Field_part* new_particle_field();
void free_particle_field(Field_part* field);

#endif
