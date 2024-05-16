#include "embed.h"
#undef EMBED
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
// #include "../double-projection2.h"
#include "../immersed.h" // IBM
#include "../interfaceforce.h" // for CD and CL
#include "curvature.h"
#include "view.h"

#define L0 15.
#define D 0.5
#define LEVEL 10

double Re;
double U0 =  1.0; // inlet velocity
double t_end = 20 [0,1];
double xi = 4.8266;
double yi = 3.0835;
coord ci = {5, 3}; // initial coordinates of cylinder
coord vc = {0, 0}; // velocity of cylinder
int j;

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

double SDF (double x, double y) {
   return - sq(x - ci.x) - sq(y - ci.y) + sq(D/2);
}

int main() {
  size(L0);
  init_grid (2 << (6));
  mu = muv;
  TOLERANCE = 1.e-7 [*]; 

  j = 10;
  Re = 1.;
  run();

  j = 11;
  Re = 2.;
  run();

  j = 12;
  Re = 5.;
  run();
  
  j = 13;
  Re = 10.;
  run();
  
  j = 14;
  Re = 20.;
  run();

  j = 15;
  Re = 40.;
  run();

  j = 20;
  Re = 1.;
  run();

  j = 21;
  Re = 2.;
  run();

  j = 22;
  Re = 5.;
  run();
  
  j = 23;
  Re = 10.;
  run();
  
  j = 24;
  Re = 20.;
  run();

  j = 25;
  Re = 40.;
  run();
}


scalar ref[];
face vector rf[];

event moving_cylinder (i++) {
  solid (airfoil, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
  solid (ref, rf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));

  if (j >= 20) {
  vector nv[];
  normal_vector(airfoil, nv);
    foreach() {
      double lambda = fabs(nv.x[]) + fabs(nv.y[]);
      if (lambda != 0.) {
        double eta = 0.065*(1 - sq(lambda)) + 0.39;
        double num = - SDF (x,y); 
        double den = lambda * eta * sqrt(2) * Delta;
        airfoil[] = den != 0? 0.5*(1 - tanh (num/den)): airfoil[];
      }
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
  coord Fp, Fmu = {0,0};
  interface_force (ref, p, u, mu, &Fp, &Fmu);
  // immersed_force (ref, &Fp);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));
 
  double E = 0;
  double E_p = 0;
  boundary ({u.x, u.y});
  scalar omega[];
  vorticity (u , omega);
  foreach(){
    double vort = omega[];
    double area = dv();
    E_p += sq(vort);
    if (airfoil[] < 1. && airfoil[] > 0){
      coord b, n;
      area *= embed_geometryo (point, &b, &n);
      vort = embed_vorticityo (point, u, b, n);
    }
    E += area*sq(vort);
  }

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, E, E_p);
	   
}


event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);
  solid (ref, rf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));

  char name[80];
  sprintf (name, "vort-%d", j);
  FILE * fp1 = fopen (name, "w");
  view (fov = 2, tx = -0.375, ty = -0.20,
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("ref", "rf", filled = 1, lw = 5);
  save (fp = fp1);
}

event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-4,3e-4,3e-4},
		 maxlevel = LEVEL, minlevel = 2);
}

event profile (t = t_end) {
  int k = 0;
  double delta = 15/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx1-%d", j); // x = 4.8125
  FILE * fv = fopen(name, "w");
  for(double i = 0; i <= 6; i += delta) {
    foreach_point (4.8125, i) {
      if (airfoil[] > 0 && airfoil[] < 1)
        k = 2.;
      else if (airfoil[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
    }
  }
  fflush (fv);
  fclose (fv);

  sprintf (name, "vprofx2-%d", j); // x = 5
  FILE * fv1 = fopen(name, "w");
  for(double i = 0; i <= 6; i += delta) {
    foreach_point (5, i) {
      if (airfoil[] > 0 && airfoil[] < 1)
        k = 2.;
      else if (airfoil[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv1, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
    }
  }
  fflush (fv1);
  fclose (fv1);

  sprintf (name, "vprofx3-%d", j); // x = 10
  FILE * fv2 = fopen(name, "w");
  for(double i = 0; i <= 6; i += delta) {
    foreach_point (10, i) {
      if (airfoil[] > 0 && airfoil[] < 1)
        k = 2.;
      else if (airfoil[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv2, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
    }
  }
  fflush (fv2);
  fclose (fv2);

  sprintf (name, "vprofy1-%d", j); // y = 3
  FILE * fv3 = fopen(name, "w");
  for(double i = 0; i <= 15; i += delta) {
    foreach_point (i, 3) {
      if (airfoil[] > 0 && airfoil[] < 1)
        k = 2.;
      else if (airfoil[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv3, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
    }
  }
  fflush (fv3);
  fclose (fv3);

  sprintf (name, "vprofy2-%d", j); // y = 3.1875
  FILE * fv4 = fopen(name, "w");
  for(double i = 0; i <= 15; i += delta) {
    foreach_point (i, 3.1875) {
      if (airfoil[] > 0 && airfoil[] < 1)
        k = 2.;
      else if (airfoil[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv4, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
    }
  }
  fflush (fv4);
  fclose (fv4);
 
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
