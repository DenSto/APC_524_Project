#ifndef PARTICLE_HPP
#define PARTICLE_HPP

typedef struct Field_part {
  double e1, e2, e3; /* Electric field components */
  double b1, b2, b3; /* Magnetic field components */
} Field_part;

typedef struct Particle {
  double x[3];  /* coordinate in X,Y,Z */	
  double v[3];  /* velocity in X,Y,Z */	

  double gamma; /* relativistic parameter gamma */

  int type;	 /* particle type */

  long my_id; 		  /* particle id */
  int initRank; 	/* initial MPI rank */

  double q;		  /* charge */
  double m;		  /* mass */
  short isGhost;    /* is particle in ghost cell? */
  short isTestParticle; /* is this particle a test particle? */

  double xo[3];  /* LAST coordinate in X,Y,Z */
  double vo[3];  /* LAST velocity in X,Y,Z */
  double dx[3]; /*last step-lengths*/

  Field_part field; /* interpolated field */

} Particle;


Particle new_particle();

Field_part new_particle_field();

#endif
