extern face vector sf; // surface solid volume fraction (0 = empty, 1 = full)
extern coord vc; // solid velocity
extern scalar airfoil;
extern scalar airfoilsf2;
face vector bf[];

face vector aF[]; // body force acceleration term


event acceleration (i++) {
  output (i, t, 10, dt);
  foreach_face()
    aF.x[] = sf.x[]*(vc.x-face_value(u.x, 0))/(dt);
  a = aF;
  output (i, t, 11, dt);
}

event end_timestep (i++) {
  output (i, t, 22, dt);
  //correction (-dt);
  foreach_face()
    aF.x[] = sf.x[]*(vc.x-face_value(u.x,0))/(dt);
  a = aF;
  //correction (dt);
  output (i, t, 23, dt);
}
/*
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



*/
