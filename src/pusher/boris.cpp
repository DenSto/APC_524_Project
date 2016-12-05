#include "boris.hpp"
#include <stdio.h>

Boris::Boris(){

}

Boris::~Boris(){
}

int Boris::Step(Particle* part, Field_part* field, double dt){
	double x,y,z,vx,vy,vz;
	double v_mx,v_my,v_mz,v_px,v_py,v_pz;
	double v_ix,v_iy,v_iz;
	double bx,by,bz,ex,ey,ez,bsqr;
	double tx,ty,tz,sx,sy,sz,tsqr;
	double x_new, y_new,z_new, vx_new, vy_new,vz_new;
	double m,q,q_p;

	x = part->x1;
	y = part->x2;
	z = part->x3;
	vx = part->v1;
	vy = part->v2;
	vz = part->v3;

	q = part->q;
	m = part->m;
	q_p = (0.5 * dt * q)/m;

	ex = field->e1;
	ey = field->e2;
	ez = field->e3;

	bx = field->b1;
	by = field->b2;
	bz = field->b3;

	tx = q_p*bx;
	ty = q_p*by;
	tz = q_p*bz;


	bsqr = bx*bx + by*by + bz*bz;
	tsqr = 2.0/(1.0 + q_p*q_p*bsqr);

	sx = tsqr*tx;
	sy = tsqr*ty;
	sz = tsqr*tz;

	v_mx = vx + q_p*ex; 
	v_my = vy + q_p*ey; 
	v_mz = vz + q_p*ez; 

	v_ix = v_mx + v_my*tz - v_mz*ty; 
	v_iy = v_my + v_mz*tx - v_mx*tz; 
	v_iz = v_mz + v_mx*ty - v_my*tx; 

	v_px = v_mx + v_iy*sz - v_iz*sy; 
	v_py = v_my + v_iz*sx - v_ix*sz; 
	v_pz = v_mz + v_ix*sy - v_iy*sx; 

	vx_new = v_px + q_p*ex; 
	vy_new = v_py + q_p*ey; 
	vz_new = v_pz + q_p*ez; 

	x_new = x + dt*vx_new;
	y_new = y + dt*vy_new;
	z_new = z + dt*vz_new;

	part->x1 = x_new;
	part->x2 = y_new;
	part->x3 = z_new;

	part->v1 = vx_new;
	part->v2 = vy_new;
	part->v3 = vz_new;

	return 1;
}

