#include "embed.h"
#include "navier-stokes/centered.h"
#include "immersed.h" // IBM
#include "interfaceforce.h" // for CD and CL
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
double t_end = 15 [0,1];

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


int main() {
  size(L0);
  init_grid (2 << (LEVEL-4));
  mu = muv;

  j = 0;
  Re = 40.;
  run();
  /*
  j++;
  Re = 20.;
  run();

  j++;
  Re = 40.;
  run();
  */
}


event init (t = 0) {
  mask(y > 6 ? top: y < -6 ? bottom : none);
}


event moving_cylinder (i++) {
  solid (airfoil, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
  boundary ((scalar *) {muv});
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
  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, E);
}


event movie (t += 0.01; t <= t_end)
{

  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = 0.5 - airfoil[];
  
  
  char name1[80];
  sprintf (name1, "%d-tracers.mp4", j);
  FILE * fp1 = fopen(name1, "w");
  output_ppm (f, fp1, n = 500, box = {{0,0},{15,6}},
	      min = -5, max = 5, linear = false, mask = m, map = cool_warm);
  fclose (fp1);

  char name2[80];
  sprintf (name2, "%d-vort.mp4", j);
  FILE * fp2 = fopen(name2, "w");
  output_ppm (omega, fp2, n = 500, box = {{0,0},{15,6}},
	      min = -5, max = 5, linear = true, mask = m, map = cool_warm); 
  fclose (fp2);
}


event snapshot (t += 2.5, t <= t_end)
{
  
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
  
  scalar omega[], m[];
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
}


event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-4,3e-4,3e-4},
		 maxlevel = LEVEL, minlevel = 2);
}


event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
