#ifndef FIELD_HPP
#define FIELD_HPP

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

Particle* new_particle();


#endif