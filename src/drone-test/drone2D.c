// #include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "utils.h"
#include "distance.h"
#include "view.h"
#if dimension == 3
  #include "lambda2.h"
#endif

scalar f[], * tracers = {f};

int LEVEL = 7;
double D = 4.6; // prop diameter
double h = 0.175; // prop height
double Re = 600000;
double U0 = 1.;
// double theta, r;
double omega = 10. [0, -1];// prop rotational speed
face vector muv[];

u.n[right]   = dirichlet(-omega*y);
u.n[left]    = dirichlet(-omega*y);
u.n[top]     = dirichlet( omega*x);
u.n[bottom]  = dirichlet( omega*x);

u.t[right]   = dirichlet( omega*fabs(x));
u.t[left]    = dirichlet(-omega*fabs(x));
u.t[top]     = dirichlet(-omega*fabs(y));
u.t[bottom]  = dirichlet( omega*fabs(y));


uf.n[right]   = dirichlet(-omega*y);
uf.n[left]    = dirichlet(-omega*y);
uf.n[top]     = dirichlet( omega*x);
uf.n[bottom]  = dirichlet( omega*x);

uf.t[right]   = dirichlet( omega*fabs(x));
uf.t[left]    = dirichlet(-omega*fabs(x));
uf.t[top]     = dirichlet(-omega*fabs(y));
uf.t[bottom]  = dirichlet( omega*fabs(y));


u.n[embed]    = dirichlet(0);
u.t[embed]    = dirichlet(0);

p[right]      = neumann (0);
p[left]       = neumann (0);
p[top]        = neumann (0);
p[bottom]     = neumann (0);

pf[right]     = neumann (0);
pf[left]      = neumann (0);
pf[top]       = neumann (0);
pf[bottom]    = neumann (0);

f[left]       = dirichlet(y > 0);

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
 refine (sq(x/0.25) + sq(y/2) < 1 && level < LEVEL);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = sq(x/0.25) + sq(y/2) - 1.;
  boundary ({phi});
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
*/
  foreach() {
    double theta = atan(y/x);
    double r = sqrt(sq(x) + sq(y));
    u.x[] = x > 0? -r*omega*sin(theta): r*omega*sin(theta);
    u.y[] = x > 0? r*omega*cos(theta): -r*omega*cos(theta);
  //  u.x[] = cs[]*u.x[];
  //  u.y[] = cs[]*u.y[];
  }

  boundary ({u});
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(uf.x[])*(D)/(Re);
  boundary ((scalar *) {muv});
}

face vector av[];

/*
event acceleration (i++) {
  coord cor = {2.0*omega,-2.0*omega};
  foreach_face(x)
    av.x[] = (cor.x*uf.y[]) + sq(omega)*x;

  foreach_face(y)
    av.y[] = cor.y*uf.x[] + sq(omega)*y;

  a = av;
}
*/

event logfile (i++; t < 5) {
 
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));

  double max_x = -1e100, max_y = -1e100, max_z = -1e100;
  foreach (reduction(max:max_x) reduction(max:max_y) reduction(max:max_z))
    foreach_dimension()
      if (fabs(u.x[]) > max_x)
        max_x = fabs(u.x[]);

  fprintf (stderr, "%d %g %g %g %g %g\n", i, t, dt, max_x, max_y, max_z);
}

/*
event graphics (i = 100)
{
  scalar omega[];
  vorticity (u, omega);
  view (fov = 20, width = 1500, height = 1500);
  squares (color = "omega", spread = 0.8, linear = true, map = cool_warm);
  vectors (u = "u", scale = 0.002);
  cells ();
  save ("fields1.png");
}
*/
/*
event movie (t += 0.01; t <= 7){
  scalar m[];
  foreach()
    m[] = cs[] - 0.5;
  static FILE * fp = fopen ("f.ppm", "w");
  output_ppm (f, fp, n = 480, linear = false,
		  min = 0, max = 0, mask = m, map = cool_warm);
}
*/

event adapt (i++) {
  double uemax = 0.5;
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(uemax),(uemax),(uemax)},
		  maxlevel = LEVEL, minlevel = (1));

  unrefine (sq(x) + sq(y) > sq(D/2) && level > 4);

}
