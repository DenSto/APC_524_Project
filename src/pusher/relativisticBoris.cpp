#include "pusher.hpp"
#include <math.h>

//! Relativistic Boris pusher
/*!
 * 	Uses the pusher described in "Simulation of beams or
 * 	plasmas crossing at relativistic velocity"
 *
 * J.-L. Vay, Phys. Plasmas 15 (5) 2007
 *
 */


class Relativistic_Boris : public Pusher {
	public:
		Relativistic_Boris();
		~Relativistic_Boris();
		int Step(Particle *part, Field_part *field, double dt);
};

Relativistic_Boris::Relativistic_Boris(){
}

Relativistic_Boris::~Relativistic_Boris(){
}

int Relativistic_Boris::Step(Particle* part, Field_part* field, double dt){
	double x,y,z,vx,vy,vz;
	double bx,by,bz,ex,ey,ez;
	double tx,ty,tz;
	double x_new, y_new,z_new, vx_new, vy_new,vz_new;
	double m,q,q_p;
	double upx,upy,upz;
	double taux,tauy,tauz,tau_sq;
	double s,sigma,ustar,udott;
	double uint_x,uint_y,uint_z;
	double uhx,uhy,uhz;
	double vXbx,vXby,vXbz;
	double upXt_x,upXt_y,upXt_z;
	double gamma,gammap_sq,gamma_int;

	x = part->x[0];
	y = part->x[1];
	z = part->x[2];
	vx = part->v[0];
	vy = part->v[1];
	vz = part->v[2];
	gamma = part->gamma;

	q = part->q;
	m = part->m;
	q_p = (0.5 * dt * q)/m;

	ex = field->e1;
	ey = field->e2;
	ez = field->e3;

	bx = field->b1;
	by = field->b2;
	bz = field->b3;

	vXbx = vy*bz - vz*by;
	vXby = vz*bx - vx*bz;
	vXbz = vx*by - vy*bx;

	uhx = gamma*vx + q_p*(ex + vXbx);
	uhy = gamma*vy + q_p*(ey + vXby);
	uhz = gamma*vz + q_p*(ez + vXbz);

	upx = uhx + q_p*ex;
	upy = uhy + q_p*ey;
	upz = uhz + q_p*ez;

	gammap_sq = 1.0 + upx*upx + upy*upy + upz*upz;

	taux = q_p*bx;
	tauy = q_p*by;
	tauz = q_p*bz;

	tau_sq = taux*taux + tauy*tauy + tauz*tauz;

	sigma = gammap_sq - tau_sq;

	ustar = upx*taux + upy*tauy + upz*tauz;

	gamma_int = sqrt(sigma*sigma + 4*(tau_sq + ustar*ustar));
	gamma_int = sqrt(0.5*(sigma + gamma_int));


	tx = taux/gamma_int;
	ty = tauy/gamma_int;
	tz = tauz/gamma_int;

	s = 1.0/(1.0 + tx*tx + ty*ty + tz*tz);
 
	udott = ustar/gamma_int;

	upXt_x = upy * tz - upz * ty;
	upXt_y = upz * tx - upx * tz;
	upXt_z = upx * ty - upy * tx;

	uint_x = s*(upx + udott*tx + upXt_x);
	uint_y = s*(upy + udott*ty + upXt_y);
	uint_z = s*(upz + udott*tz + upXt_z);

	vx_new = uint_x/gamma_int; 
	vy_new = uint_y/gamma_int; 
	vz_new = uint_z/gamma_int; 

	x_new = x + dt*vx_new;
	y_new = y + dt*vy_new;
	z_new = z + dt*vz_new;

	//Update last particle position.
	part->xo[0] = x;
	part->xo[1] = y;
	part->xo[2] = z;

	//Update particle step length
	part->dx[0] = x_new-x;
	part->dx[1] = y_new-y;
	part->dx[2] = z_new-z;


	//Update last particle velocity.
	part->vo[0] = vx;
	part->vo[1] = vy;
	part->vo[2] = vz;

	part->x[0] = x_new;
	part->x[1] = y_new;
	part->x[2] = z_new;

	part->v[0] = vx_new;
	part->v[1] = vy_new;
	part->v[2] = vz_new;

	part->gamma = gamma_int;

	return 1;
}

