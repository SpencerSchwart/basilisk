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
#line 1 "movingairfoil.c"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "../immersed.h" // IBM
#include "../interfaceforce.h" // for CD and CL
#include "curvature.h"
#include "view.h"

#define MAX_THICKNESS 0.12
#define CHORD_LENGTH 1
#define L0 20.
#define Re (6000000.)
#define LEVEL 13

double U0 =  1.0; // inlet velocity
double rr = 1.1019*sq(MAX_THICKNESS); // Radius of leading edge
double theta_p = 0.0872665; // aoa = 5 degrees

coord vc = {0.,0.}; // the velocity of the cylinder
coord ci = {5, 10}; // initial coordinates of airfoil
coord cr = {0.25*(CHORD_LENGTH), 0.}; // center of rotation

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

u.n[top]   = neumann (0);
u.n[bottom] = neumann (0);

#define airfoil_thickness(x) (5 * MAX_THICKNESS * ((0.2969*(sqrt(x)))	\
						   -(0.1260*x)		\
						   -(0.3516*(sq(x)))	\
						   +(0.2843*(cube(x)))	\
						   -(0.1036*(pow(x, 4.)))))


void airfoil_shape (scalar c, face vector f, double theta, vertex scalar phii = {0})
{
  double yt, yc = 0;
  vertex scalar phi = automatic (phii);
  
  double chord = CHORD_LENGTH;
  
  foreach_vertex() {
    
    double XX = cr.x + (x - ci.x)*cos (theta) - (y - ci.y)*sin (theta);
    double YY = cr.y + (x - ci.x)*sin (theta) + (y - ci.y)*cos (theta);

    if (XX < 0) {
      // leading edge approximation
      phi[] = - sq(XX-rr) - sq(YY) + sq(rr);
    }
    else if (XX >= 0. && XX <= chord) {
      // basic airfoil thickness
      yt = airfoil_thickness(XX);
      
      phi[] = YY > yc? - sq(YY) + sq(yt): - sq(YY) + sq(yt);
    }
    else
      phi[] = -1.;
  } 
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f);
  c.refine = c.prolongation = fraction_refine;
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


double SDF (double xx, double yy) {
  double XX = cr.x + (xx - ci.x)*cos (theta_p) - (yy - ci.y)*sin (theta_p);
  double YY = cr.y + (xx - ci.x)*sin (theta_p) + (yy - ci.y)*cos (theta_p);

  if (XX < 0) {
      // leading edge approximation
      return - sq(XX-rr) - sq(YY) + sq(rr);
    }
  else if (XX >= 0. && XX <= CHORD_LENGTH) {
      // basic airfoil thickness
      double yt = airfoil_thickness(XX);
      return YY > 0? - sq(YY) + sq(yt): - sq(YY) + sq(yt);
    }
  else
      return -1.;
}	


int main(){
  size(L0);
  init_grid (2 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5 [*];
  run();
}


event init (t = 0) {
  // initial mesh refinement
  astats ss;
  int ic = 0;
  do {
    ic++;
    airfoil_shape (airfoil, sf, theta_p);
    ss = adapt_wavelet ({airfoil}, (double[]) {1.e-30},
			maxlevel = LEVEL, minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
  airfoil_shape(airfoil, sf, theta_p);
}

scalar ref[];
face vector rf[];
vertex scalar phia[];

event moving_cylinder (i++) {
  airfoil_shape (airfoil, sf, theta_p, phia);
  airfoil_shape (ref, rf, theta_p);

  foreach()
    if (airfoil[] > 0)  {
      vector nv[];
      normal_vector(ref, nv);
      double lambda = fabs(nv.x[]) + fabs(nv.y[]);
      double eta = 0.065*(1 - sq(lambda)) + 0.39;
      double num = - SDF (x,y);
      double den = lambda * eta * sqrt(2) * Delta;
      double scale = 15;
      double vof = airfoil[];
      airfoil[] = den != 0? 0.5*(1 - tanh ((num/den)*scale)):airfoil[];
      // fprintf (stderr, "i=%d t=%g n.x=%g n.y=%g lam=%g eta=%g vof_b=%g vof_a=%g x=%g y=%g delta=%g num=%g den=%g tanh(%g)\n",
// 		      i, t, nv.x[], nv.y[], lambda, eta, vof, airfoil[], x, y, Delta, num, den, (num/den)*scale);
    }


}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(CHORD_LENGTH)/(Re);
  boundary ((scalar *) {muv});
}


event logfile (i++; t <= 15){

  coord Fp, Fmu;
  interface_force (airfoil, p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(CHORD_LENGTH));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(CHORD_LENGTH));
  fprintf (stderr, "%d %g %d %d %d %d %g %g\n",
	   i, t, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
  
}


/*
event movie (t += 1e-2; t <= 10)
{
  scalar omega[];
  vorticity (u, omega);
  view (fov = 17, camera = "front",
	tx = -0.5, ty = -0.5, bg = {1,1,1},
	width = 3200, height = 3200);
  clear();
  draw_vof ("airfoil", "sf",filled = 1, lw = 3);
  squares ("u.x", map = cool_warm);
  //squares ("omega", map = cool_warm);;
  save ("vorticity.mp4");
}
*/

#if DUMP
event snapshot (t += 1e-1; t <= 5)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}
#endif

event coordinates (t=0.1) {
  char name[80];
  char name1[80];
  sprintf (name, "out-%g", t);
  sprintf (name1, "out-facets");
  FILE * fp = fopen (name, "w");
  FILE * fp1 = fopen (name1, "w");
  foreach() {
    if (airfoil[] > 0 && airfoil[] < 1)
      fprintf (fp, "%g %g\n", x, y);
  }
  fprintf (fp, "\n");
  fclose (fp);
  output_facets (airfoil, fp1, sf);
  fclose (fp1);

}


event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-4,3e-4,3e-4},
		 maxlevel = LEVEL, minlevel = LEVEL - 6);
}

#endif
