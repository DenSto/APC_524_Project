#ifndef PARTICLE_HPP
#define PARTICLE_HPP

typedef struct Field_part {
  double e1, e2, e3; /* Electric field components */
  double b1, b2, b3; /* Magnetic field components */
} Field_part;

typedef struct Particle {
  double x[3];  /* coordinate in X,Y,Z */	
  double v[3];  /* velocity in X,Y,Z */	

  double xo[3];  /* LAST coordinate in X,Y,Z */
  double vo[3];  /* LAST velocity in X,Y,Z */

  double q;		  /* charge */
  double m;		  /* mass */

  int my_id; 		  /* particle id */
  short isGhost;    /* is particle in ghost cell? */

  Field_part field; /* interpolated field */

} Particle;


Particle new_particle();

Field_part new_particle_field();

#endif
