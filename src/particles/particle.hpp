#ifndef PARTICLE_HPP
#define PARTICLE_HPP

typedef struct Field_part {
  double e1, e2, e3; /* Electric field components */
  double b1, b2, b3; /* Magnetic field components */
} Field_part;

typedef struct Particle {
  double x1,x2,x3;  /* coordinate in X,Y,Z */	
  double v1,v2,v3;  /* velocity in X,Y,Z */	

  double xo1,xo2,xo3;  /* LAST coordinate in X,Y,Z */
  double vo1,vo2,vo3;  /* LAST velocity in X,Y,Z */

  double q;		  /* charge */
  double m;		  /* mass */

  int my_id; 		  /* particle id */
  short isGhost;    /* is particle in ghost cell? */

  Field_part* field; /* interpolated field */

} Particle;


Particle* new_particle();
void free_particle(Particle* part);

Field_part* new_particle_field();
void free_particle_field(Field_part* field);

#endif
