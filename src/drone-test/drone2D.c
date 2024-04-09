#include "embed.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "view.h"

int LEVEL = 7;
double D = 4.6; // prop diameter
double h = 0.175; // prop height
double Re = 600000;
double U0 = 1.;
// double theta, r;
double omega = 10. [0, -1];// prop rotational speed
face vector muv[];


u.n[right]   = dirichlet( omega*y);
u.n[left]    = dirichlet( omega*y);
u.n[top]     = dirichlet(-omega*x);
u.n[bottom]  = dirichlet(-omega*x);

u.t[right]   = dirichlet(-omega*fabs(y));
u.t[left]    = dirichlet( omega*fabs(y));
u.t[top]     = dirichlet( omega*fabs(x));
u.t[bottom]  = dirichlet(-omega*fabs(x));

/*
uf.n[right]   = dirichlet( omega*y);
uf.n[left]    = dirichlet( omega*y);
uf.n[top]     = dirichlet(-omega*x);
uf.n[bottom]  = dirichlet(-omega*x);

uf.t[right]   = dirichlet(-omega*x);
uf.t[left]    = dirichlet(-omega*x);
uf.t[top]     = dirichlet( omega*y);
uf.t[bottom]  = dirichlet( omega*y);
*/
u.n[embed]    = dirichlet(0);
u.t[embed]    = dirichlet(0);

/*
scalar psi[];
psi[right]  = dirichlet(-0.5*omega*(sq(x) + sq(y)));
psi[left]   = dirichlet(-0.5*omega*(sq(x) + sq(y)));
psi[top]    = dirichlet(-0.5*omega*(sq(x) + sq(y)));
psi[bottom] = dirichlet(-0.5*omega*(sq(x) + sq(y)));
*/

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
  double L0 = 5;
  size (L0);
  origin (-L0/2, -L0/2,-L0/2);
  mu = muv;
  run();
}

event init (t = 0) {
/*
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("/home/spencer/basilisk/CTFL/files/prop_test5.stl", "r");
    if (fp == NULL)
      fprintf(stderr, "STL file is NULL\n");
    fraction_from_stl (cs, fs, fp);
    fclose (fp);
  }
  */
/* 
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = fabs(x) < 0.25 && fabs(y) < 1? 0:1;	  
  }
  boundary ({phi});
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  */
   foreach()
    {
	double theta = atan(y/x);
	
	double r = sqrt(sq(x) + sq(y));
	u.x[] = x > 0? r*omega*sin(theta):-r*omega*sin(theta);
	u.y[] = x > 0? -r*omega*cos(theta): r*omega*cos(theta);
    }

}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(uf.x[])*(D)/(Re);
  boundary ((scalar *) {muv});
}

face vector av[];

event acceleration (i++) {
  coord cor = {-2.0*omega, 2.0*omega};
  foreach_face(x)
    av.x[] = (cor.x*uf.y[]) + sq(omega)*x;

  foreach_face(y)
    av.y[] = cor.y*uf.x[] + sq(omega)*y;

  a = av;
}


event logfile (i++, t <= 5) {
 /*
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
  */
 /*
  double E = 0; 
  boundary ({u.x, u.y});
  scalar omegav[];
  vorticity (u , omegav);
  foreach(){
    double vort = omegav[];
    double area = dv();
    if (cs[] < 1. && cs[] > 0){
      coord b, n;
      area *= embed_geometry (point, &b, &n);
      vort = embed_vorticity (point, u, b, n);
    }
    E += area*sq(vort);
  }
 
 */
  fprintf (stderr, "%d %g\n", i, t);
}


event graphics (i = 1)
{
  scalar omega[];
  vorticity (u, omega);
  view (fov = 20, width = 1500, height = 1500);
  squares (color = "omega", spread = 0.8, linear = true, map = cool_warm);
  vectors (u = "u", scale = 0.002);
  cells ();
  save ("fields1.png");
}


event movie (t += 0.01){
}

event adapt (i++) {
  double uemax = 0.001;
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(1.e-3),(1.e-3)},
		  maxlevel = LEVEL, minlevel = (1));
}
