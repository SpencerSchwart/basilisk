#ifndef BASILISK_HEADER_25
#define BASILISK_HEADER_25
#line 1 "./../immersed.h"
extern face vector sf; // surface solid volume fraction (0 = empty, 1 = full)
extern coord vc; // solid velocity
extern scalar airfoil;
extern scalar airfoilsf2;
face vector bf[];

face vector aF[]; // body force acceleration term


event acceleration (i++) {
  foreach_face() {
    // double ff = (airfoilsf[] + airfoilsf[-1])/2;
    aF.x[] = face_value(airfoil,0)*(vc.x-face_value(u.x, 0))/(dt);
    // bf.x[] += sf.x[]*(vc.x-face_value(u.x, 0))/(dt);
  }
  a = aF;
}

event end_timestep (i++) {
  correction (-dt);
  foreach_face()
    aF.x[] = face_value(airfoil,0)*(vc.x-face_value(u.x,0))/(dt);
  a = aF;
  correction (dt);
}

void imersed_force (scalar c, coord * F) {

  coord Fi = {0}; // intermediate force
  vector aC[]; // body force acceleration term (centered)

  foreach() {
    if (c[] > 0. && c[] < 1.) {
      coord n = interface_normal (point, c), p;
      double alphaP = plane_alpha (c[], n);
      double area = pow(Delta, dimension - 1)*plane_area_center (n, alphaP, &p);
      aC.x[] = (aF.x[] + aF.x[1])/2; // average surface velocity to get center
      Fi.x += -(aC.x[]*n.x*area)/alpha.x[];
      aC.y[] = (aF.y[0, 0] + aF.y[0, 1])/2;
      Fi.y += -(aC.y[]*n.y*area)/alpha.y[];
    }
  }   
  *F = Fi; 
}





#endif
