// #include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "view.h"

int LEVEL = 8;
double D = 4.6; // prop diameter
double h = 0.175; // prop height
double Re = 600000;
double U0 = 1.;
double theta, r;
double omega = 1. [0, -1]; // prop rotational speed
scalar cad[];
face vector muv[];

// u.n[bottom] = dirichlet(0.);
p[bottom]   = neumann(0.);
pf[bottom]  = neumann(0.);
// u.n[top]  = neumann(0.);
p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);
p[right] = neumann(0.);
pf[right] = neumann(0.);

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);
			

void fraction_from_stl (scalar cs, face vector fs, FILE * fp) {
  coord * p = input_stl (fp); // import CAD
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;

  // STL points to domain points
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){5e-4}, LEVEL, 5).nf);

  // VOF calculation
  vertex scalar phi[];
  foreach_vertex(){
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
}

int main() {
  fprintf (stderr, "main start\n");
  init_grid (32);
  // mu = muv;
  double L0 = 5;
  size (L0);
  origin (-L0/2, -L0/2,-L0/2);
  
  for (scalar s in {cad})
    s.refine = s.prolongation = fraction_refine;
  fprintf (stderr, "before BC\n");

  
  
  fprintf (stderr, "after BC\n");
  run();
}

event init (t = 0) {

  /*
  fprintf(stderr, "init start\n");
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("/home/spencer/basilisk/src/examples/prop_test5.stl", "r");
    if (fp == NULL)
      fprintf(stderr, "STL file is NULL\n");
    fraction_from_stl (cs, fs, fp);
    fclose (fp);
  }
  */

  /*
  foreach()
    foreach_dimension()
    u.x[] = cs[]*1.;
  */

  // && abs(z) <= (h/2) && cs[] < 1

  foreach()
    {
      if (abs(x) <= D/2 && abs(y) <= D/2)
	if (sq(x) + sq(y) < sq(D/2)) {
	  theta = atan(y/x);
	  r = sqrt(sq(x) + sq(y));
	  u.x[] = x > 0? -1*r*omega*sin(theta):r*omega*sin(theta);
	  u.y[] = x > 0? r*omega*cos(theta): -1*r*omega*cos(theta);
	}
    }
  
}

/*
u.t[bottom] = dirichlet(x > 0?-1*(sqrt(sq(x) + sq(y)))*omega*sin(atan(y/x)):
			(sqrt(sq(x) + sq(y)))*omega*sin(atan(y/x)));
u.n[bottom] = dirichlet(x > 0? (sqrt(sq(x) + sq(y)))*omega*cos(atan(y/x)):
			-1*(sqrt(sq(x) + sq(y)))*omega*cos(atan(y/x)));
*/


event properties (i++) {
  /*
  foreach_face()
    muv.x[] = fm.x[]*(uf.x[])*(D)/(Re);
  boundary ((scalar *) {muv});
  */

  
  foreach_boundary(bottom) {
    theta = atan(y/x);
    r = sqrt(sq(x) + sq(y));
    u.x[] = x > 0? -1*r*omega*sin(theta):r*omega*sin(theta);
    u.y[] = x > 0? r*omega*cos(theta): -1*r*omega*cos(theta);
  }
  foreach_boundary(top) {
    theta = atan(y/x);
    r = sqrt(sq(x) + sq(y));
    u.x[] = x > 0? -1*r*omega*sin(theta):r*omega*sin(theta);
    u.y[] = x > 0? r*omega*cos(theta): -1*r*omega*cos(theta);
  }
  foreach_boundary(left) {
    theta = atan(y/x);
    r = sqrt(sq(x) + sq(y));
    u.x[] = x > 0? -1*r*omega*sin(theta):r*omega*sin(theta);
    u.y[] = x > 0? r*omega*cos(theta): -1*r*omega*cos(theta);
  }
  foreach_boundary(right) {
    theta = atan(y/x);
    r = sqrt(sq(x) + sq(y));
    u.x[] = x > 0? -1*r*omega*sin(theta):r*omega*sin(theta);
    u.y[] = x > 0? r*omega*cos(theta): -1*r*omega*cos(theta);
  }

}

event logfile (i++, t <= 5) {

  /*
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
  */
  fprintf (stderr, "%d %g %g %g\n", i, t);
}


event graphics (i = 0)
{
  scalar omega[];
  vorticity (u, omega);
  view (fov = 20, width = 1500, height = 1500);
  squares (color = "omega", spread = 0.8, linear = true, map = cool_warm);
  vectors (u = "u", scale = 0.01);
  cells ();
  save ("fields0.png");
}


event adapt (i++) {
  double uemax = 0.001;
  adapt_wavelet ({cs,u},
		 (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
}



