#include "boris.hpp"
#include <stdio.h>
#include "../globals.hpp"

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

  x = part->x[0];
  y = part->x[1];
  z = part->x[2];
  vx = part->v[0];
  vy = part->v[1];
  vz = part->v[2];

  q = part->q;
  if(debug>3)fprintf(stderr,"Charge inside = %f\n",q);

  m = part->m;
  if(debug>3)fprintf(stderr,"Mass inside = %f\n",m);

  q_p = (0.5 * dt * q)/m; 
  q_p *= UNIT_ACC; // multiply unit of acceleration
  if(debug>3)fprintf(stderr,"q_p inside = %f\n",q_p);

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

  if(debug>3){
      fprintf(stderr,"*************************\n");
      fprintf(stderr,"old: x=%f,y=%f,x=%f,vx=%f,vy=%f,vz=%f\n",
               x,y,z,vx,vy,vz);
      fprintf(stderr,"new: x=%f,y=%f,x=%f,vx=%f,vy=%f,vz=%f\n",
               x_new,y_new,z_new,vx_new,vy_new,vz_new);
      fprintf(stderr,"*************************\n");
  }

  return 1;
}

