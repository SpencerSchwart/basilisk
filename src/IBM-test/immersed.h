#ifndef BASILISK_HEADER_23
#define BASILISK_HEADER_23
#line 1 "./../immersed.h"
extern coord vc; // solid velocity
extern scalar airfoil;
extern face vector sf;
vector aF[]; // body force acceleration term

/*
event acceleration (i++) {
  foreach_face()
    aF.x[] = face_value(airfoil,0)*(vc.x-face_value(u.x,0))/(dt);
  a = aF;
}
*/

event end_timestep (i++) {

  foreach()
    foreach_dimension()
      aF.x[] = airfoil[]*(vc.x-u.x[])/(dt);

  correction(-dt);
  foreach()
    foreach_dimension()
      g.x[] += fm.x[]*aF.x[];
  correction(dt);

}







void normal_vector (scalar s, vector nv) {
  foreach()
    foreach_dimension()
      nv.x[] = (s[1] - s[-1])/(2.*Delta);

  foreach() {
    double mag = 0;
    foreach_dimension()
      mag += sq(nv.x[]);
    mag = sqrt(mag);
    if (mag > 0)
      foreach_dimension()
	nv.x[] /= mag;
  }
}

double int_area (Point point, coord * b, coord * n) {
  *n = facet_normal (point, airfoil, sf);
  double alpha = plane_alpha (airfoil[], *n);
  double area = plane_area_center (*n, alpha, b);
  normalize (n);
  return area;
}


void immersed_force (scalar c, coord * F) {
  coord Fi = {0, 0}; // intermediate force
  vector nv[];
  // normal_vector (c, nv);
  foreach() {
    if (c[] > 0. && c[] < 1.) {
      coord n,b;
      double area = int_area (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      double Fn_x = aF.x[]*n.x + aF.y[]*n.x;
      double Fn_y = aF.x[]*n.y + aF.y[]*n.y;
      foreach_dimension()
        Fi.x += (Fn_x*area);
    }
  }
  *F = Fi; 
}



#endif
