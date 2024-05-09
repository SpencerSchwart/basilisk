#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "utils.h"
#include "distance.h"
#include "lambda2.h"
#include "view.h"

int LEVEL = 8;
double D = 1.05; // prop diameter
double h = 0.05; // prop height
double Re = 600000;
double omega = 20. [0, -1];// prop rotational speed
face vector muv[];
int j;

u.n[right]    = dirichlet( omega*y);
u.t[right]    = dirichlet(-omega*fabs(x));
uf.n[right]   = dirichlet( omega*y);
uf.t[right]   = dirichlet(-omega*fabs(x));
p[right]      = neumann (0);
pf[right]     = neumann (0);
u.r[right]    = dirichlet(0);

u.n[left]     = dirichlet( omega*y);
u.t[left]     = dirichlet( omega*fabs(x));
uf.n[left]    = dirichlet( omega*y);
uf.t[left]    = dirichlet( omega*fabs(x));
p[left]       = neumann (0);
pf[left]      = neumann (0);
u.r[left]     = dirichlet(0);

u.n[top]      = dirichlet(-omega*x);
u.t[top]      = dirichlet(0);
uf.n[top]     = dirichlet(-omega*x);
//uf.t[top]     = dirichlet( omega*fabs(y));
p[top]        = neumann (0);
pf[top]       = neumann (0);
u.n[top]      = dirichlet( omega*fabs(y));

u.n[bottom]   = dirichlet(-omega*x);
u.r[bottom]   = dirichlet(0);
uf.n[bottom]  = dirichlet(-omega*x);
//uf.t[bottom]  = dirichlet(-omega*fabs(y));
p[bottom]     = neumann (0);
pf[bottom]    = neumann (0);
u.t[bottom]   = dirichlet(-omega*fabs(y));

// Ceiling
u.n[front]    = dirichlet(0);
u.t[front]    = dirichlet( omega*y);
u.r[front]    = dirichlet(-omega*x);

// Outlet
u.n[back]     = neumann(0);
u.t[back]     = dirichlet( omega*y);
p[back]       = neumann (0);
pf[back]      = neumann (0);
u.r[back]     = dirichlet(-omega*x);

u.n[embed]    = dirichlet(0);
u.t[embed]    = dirichlet(0);


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
  foreach_vertex(){;
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
}


int main() {
  fprintf (stderr, "main start\n");
  init_grid (32);
  double L0 = 2, d;
  size (L0);
  mu = muv;

  j = 0;
  d = 0.9;
  origin (-L0/2, -L0/2, -(L0-(D*d)));
  run();

  j++;
  d = 0.1;
  origin (-L0/2, -L0/2, -(L0-(D*d)));
  run();
 
  j++;
  d = 0.2;
  origin (-L0/2, -L0/2, -(L0-(D*d)));
  run();

  j++;
  d = 0.4;
  origin (-L0/2, -L0/2, -(L0-(D*d)));
  run();

  j++;
  d = 0.6;
  origin (-L0/2, -L0/2, -(L0-(D*d)));
  run();

  j++;
  d = 0.15;
  origin (-L0/2, -L0/2, -(L0-(D*d)));
  run();
}


event init (t = 0) {
 
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("/home/spencer/basilisk/CTFL/files/prop_test1.stl", "r");
    if (fp == NULL)
      fprintf(stderr, "STL file is NULL\n");
    else {
      fprintf(stderr, "STL is NOT null\n");
      fraction_from_stl (cs, fs, fp);
      foreach()
        cs[] = 1 - cs[];
      foreach_face()
	fs.x[] = 1 - fs.x[];
      fclose (fp);
    }
  }

  foreach() {
    u.x[] = omega * y; 
    u.y[] = -omega * x;
    u.z[] = 0.;
  }
  boundary ({u});
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(uf.x[])*(D)/(Re);
  boundary ((scalar *) {muv});
}


event logfile (i++; t <= 10) {

  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x);
  double CLy = (Fp.y + Fmu.y);
  double CLz = (Fp.z + Fmu.z);

  double max_x = -1e100, max_y = -1e100, max_z = -1e100;
  foreach (reduction(max:max_x) reduction(max:max_y) reduction(max:max_z))
    foreach_dimension()
      if (fabs(u.x[]) > max_x)
        max_x = fabs(u.x[]);

  fprintf (stderr, "%d %g %d %g %g %g %g %g %g %g\n", i, t, j, dt, max_x, max_y, max_z, CD, CLy, CLz);
}

scalar l2[];
event movie (t += 0.01) {
  view (fov = 30, quat = {0.515965,0.140691,0.245247,0.808605}, width = 1024, height = 768);

  clear();
  draw_vof ("cs");
  lambda2 (u, l2);
  isosurface ("l2", 1);
  if (j == 0)
    save ("0-movie.mp4");
  if (j == 1)
    save ("1-movie.mp4");
  if (j == 2)
    save ("2-movie.mp4");
  if (j == 3)
    save ("3-movie.mp4");
  if (j == 4)
    save ("4-movie.mp4");
  if (j == 5)
    save ("5-movie.mp4");
}

event snapshot (t += 1) {
  char name[80];
  sprintf (name, "%d-dump-%g", j, t);
  dump (name); 
}

event image (t += 1.25) {
  char name[80];
  sprintf (name, "%d-overview-%g.png", j, t);
  view (fov = 30, quat = {0.515965,0.140691,0.245247,0.808605}, width = 1024, height = 768);
  clear();
  draw_vof ("cs", "fs");
  lambda2 (u, l2);
  isosurface ("l2", 1);
  save (name);

  view (fov = 20, tx = -0.45, camera = "right", width = 1024, height = 1024);
  sprintf (name, "%d-vorticity-%g.png", j, t);
  clear();
  draw_vof ("cs", "fs"); 
  scalar omega[];
  vorticity (u, omega);
  squares ("omega", n = {1,0,0}, map = cool_warm);
  save (name);

  sprintf (name, "%d-pressure-%g.png", j, t);
  clear();
  draw_vof ("cs", "fs"); 
  squares ("p", n = {1,0,0}, map = cool_warm);
  save (name);

  sprintf (name, "%d-yzcells-%g.png", j, t);
  clear();
  draw_vof ("cs", "fs"); 
  squares ("omega", n = {1,0,0}, map = cool_warm);
  cells (n = {1,0,0});
  save (name);

  view (fov = 20, camera = "front");
  sprintf (name, "%d-xycells-%g.png", j, t);
  clear();
  draw_vof ("cs", "fs"); 
  squares ("omega", n = {0,0,1}, map = cool_warm);
  cells (n = {0,0,1});
  save (name);

  sprintf (name, "%d-zvelocity-%g.png", j, t);
  clear();
  draw_vof ("cs", "fs"); 
  squares ("u.z", n = {0,0,1}, map = cool_warm);
  cells (n = {0,0,1});
  save (name);

  sprintf (name, "%d-outlet-%g.png", j, t);
  clear();
  draw_vof ("cs", "fs"); 
  squares ("u.z", n = {0,0,1}, alpha = -2.3, map = cool_warm);
  cells (n = {0,0,1}, alpha = -2.3);
  save (name); 
}

event adapt (i++) {
  double uemax = 0.1;
  adapt_wavelet ({cs,u}, (double[]) {1.e-2, (uemax),(uemax),(uemax)},
		  maxlevel = LEVEL, minlevel = (2));
  unrefine (sq(x) + sq(y) > sq(D/2) && level > 4);
  // unrefine (z < -2.00 && level > 3);
}
