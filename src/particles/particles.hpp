#ifndef FIELD_HPP
#define FIELD_HPP

typedef struct Particle {
	double x1,x2,x3;  /* coordinate in X,Y,Z */	
	double v1,v2,v3;  /* velocity in X,Y,Z */	
	
	int my_id; 		  /* particle id */
	short isGhost;    /* is particle in ghost cell? */

} Particle;

#endif
