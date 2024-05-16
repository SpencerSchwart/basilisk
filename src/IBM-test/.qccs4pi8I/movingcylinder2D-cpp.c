@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "movingcylinder2D.c"
#include "embed.h"
#undef EMBED
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "../immersed.h" // IBM
#include "../interfaceforce.h" // for CD and CL
#include "curvature.h"
#include "view.h"
#include "tracer.h"

#define L0 15.
#define D 0.5
double Re;
#define LEVEL 10

double U0 =  1.0; // inlet velocity
coord ci = {5, 3}; // initial coordinates of airfoil
coord vc = {0, 0};
int j;
double t_end = 20 [0,1];

double xi = 4.8266;
double yi = 3.0835;

scalar f[];
scalar * tracers = {f};
scalar airfoil[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet ((U0));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);
p[top] = neumann (0);
u.n[bottom] = neumann (0);
p[bottom] = neumann (0);

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


int main() {
  size(L0);
  init_grid (2 << (6));
  mu = muv;
  TOLERANCE = 1.e-5 [*]; 
  
  j = 0;
  Re = 40.;
  run();
/*
  j++;
  Re = 2.;
  run();

  j++;
  Re = 5.;
  run();
  
  j++;
  Re = 10.;
  run();
  
  j++;
  Re = 20.;
  run();

  j++;
  Re = 40.;
  run();
  */
}

scalar reference[];
face vector rf[];
event moving_cylinder (i++) {
  solid (airfoil, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
  solid (reference, rf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
  /*
  foreach()
    foreach_dimension()
      u.x[] = airfoil[]*vc.x + (1. - airfoil[])*u.x[];
  boundary ((scalar *){u});
  */
 
  foreach() {
    if (airfoil[] > 0 && airfoil[] < 1) {
      vector nv[];
      normal_vector(reference, nv);
      double lambda = fabs(nv.x[]) + fabs(nv.y[]);
      double eta = 0.065*(1 - sq(lambda)) + 0.39;
      double k = 1;
      double vof = airfoil[];
      double num = sqrt(sq(x - ci.x) + sq(y - ci.y)) - (D/2);
      //fprintf (stderr, "alpha i=%d t=%g n.x=%g n.y=%g lam=%g eta=%g vof=%g x=%g y=%g delta=%g num=%g\n",i,t,
      //         nv.x[], nv.y[], lambda, eta, vof, x, y, Delta, num);
     //airfoil[] = 0.5 * (1 - tanh (reference[] / (lambda * eta * sqrt(2) * Delta)));
     airfoil[] = 0.5*(1 - tanh ( num / (lambda * eta * sqrt(2) * Delta)));
      k++;
      vof = airfoil[];
      //fprintf (stderr, "alpha k=%g i=%d t=%g n.x=%g n.y=%g lam=%g eta=%g vof=%g x=%g y=%g delta=%g num=%g\n",k,i,t,
      //       nv.x[], nv.y[], lambda, eta, vof, x, y, Delta, num);
    }
  }
}

event init (t = 0) {
  mask(y > 6 ? top: y < -6 ? bottom : none);
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
 // boundary ((scalar *) {muv});
}


event logfile (i++){

  coord Fp, Fmu;
  interface_force (airfoil, p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));
 
  double E = 0;
  boundary ({u.x, u.y});
  scalar omega[];
  vorticity (u , omega);
  foreach(){
    double vort = omega[];
    double area = dv();
    if (airfoil[] < 1. && airfoil[] > 0){
      coord b, n;
      area *= embed_geometryo (point, &b, &n);
      vort = embed_vorticityo (point, u, b, n);
    }
    E += area*sq(vort);
  }
  int counter = 0;
  double u_avg_x = 0., u_avg_y = 0., p_avg = 0.;
  foreach()
    if (airfoil[] == 1) {
      p_avg += p[];
      foreach_dimension()
        u_avg_x += u.x[];
       counter++;	
    }
  u_avg_x /= counter;
  u_avg_y /= counter;
  p_avg /= counter;

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, u_avg_x, u_avg_y, p_avg, CD, CL, E);
}


event movie (t += 0.05; t <= t_end)
{
/*
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = 0.5 - airfoil[];
  

  char name[80];

  sprintf (name, "%d-tracers.mp4", j);
  FILE * fp1 = fopen(name, "w");
  output_ppm (f, fp1, n = 500, box = {{0,0},{15,6}},
	      min = -5, max = 5, linear = false, mask = m, map = cool_warm);
  fclose (fp1);

  sprintf (name, "%d-vort.mp4", j);
  static FILE * fp2 = fopen(name, "w");
  output_ppm (omega, fp2, n = 500, box = {{4,2},{15,4}},
	      min = -5 , max = 5, linear = true, mask = m, map = cool_warm); 
  // fclose (fp2);

  view (fov = 0.75, tx = -0.3325, ty = -0.20,
		  width = 800, height = 600);
  clear();
  draw_vof ("airfoil", "sf", filled = -1, lw = 5);
  vectors("u", scale = 0.02);
  squares("u.x", min = -0.05, max = 0.05, map = cool_warm);
  save ("vinside.mp4");


  clear();
  draw_vof ("airfoil", "sf", filled = -1, lw = 5);
  vectors("u", scale = 0.02);
  squares("p", map = cool_warm);
  save ("0-pinside.mp4");
  */
}


scalar psi[];
event snapshot (t += 5, t <= t_end)
{

  scalar omega[];
  vertex scalar stream[];

  vorticity (u, omega);
  // stream lines
  psi[bottom] = dirichlet (0);
  psi[top] = dirichlet (0);
  psi[left] = dirichlet (0);
  psi[right] = dirichlet (0);
  poisson (psi, omega);
  boundary ({psi});
  foreach_vertex()
    stream[] = (psi[0,-1] + psi[-1,-1]
		+ psi[] + psi[-1])/4;
  /*
  char name[80];
  sprintf (name, "%d-dump-%d", j, i);
  dump (file = name);
  
  scalar m[];
  vorticity (u, omega);
  foreach()
    m[] = 0.5 - airfoil[];
  
  char name1[80];
  sprintf (name1, "%d-vort-%g.png", j, t);
  FILE * fp1 = fopen(name1, "w");
  output_ppm (omega, fp1, n = 1000, box = {{0,0},{15,6}},
	      min = -5, max = 5, linear = true, mask = m, map = cool_warm); 
  fclose (fp1);
  
  char name2[80];
  sprintf (name2, "%d-velocity-%g.png", j, t);
  FILE * fp2 = fopen(name2, "w");
  output_ppm (u.x, fp2, n = 1000, box = {{0,0},{15,6}},
	      linear = true, mask = m, map = cool_warm); 
  fclose (fp2);

  char name3[80];
  sprintf (name3, "%d-pressure-%g.png", j, t);
  FILE * fp3 = fopen(name3, "w");
  output_ppm (p, fp3, n = 1000, box = {{0,0},{15,6}},
	      linear = true, mask = m, map = cool_warm); 
  fclose (fp3);
  */

  char name[80];
  sprintf (name, "%d-stream-lines-%g", j, t);
  FILE * fp4 = fopen (name, "w");
  view (fov = 2, tx = -0.375, ty = -0.20,
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("airfoil", "sf", filled = 1, lw = 5);
  save (fp = fp4);

  /*
  sprintf (name, "%d-vinside-%g", j, t);
  FILE * fp5 = fopen (name, "w");
  view (fov = 0.75, tx = -0.3325, ty = -0.20,
		  width = 800, height = 600);
  clear();
  draw_vof ("airfoil", "sf", filled = -1, lw = 5);
  vectors("u", scale = 0.02);
  squares("u.x", min = -0.05, max = 0.05, map = cool_warm);
  save (fp = fp5);
  */
}


event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-4,3e-4,3e-4},
		 maxlevel = LEVEL, minlevel = 2);
}

event profile (t = 0.1) {
  static FILE * fp = fopen("vprof", "w");
  int count = 0;
  foreach() {
    double xx = x;
    foreach_point (xx, 5.) {
    count += 1;
    }
  }
  fprintf (fp, "%d", count);
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

#endif
