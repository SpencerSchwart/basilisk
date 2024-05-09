#include "embed.h"
#undef EMBED
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

u.n[top] = neumann (0);
p[top] = dirichlet (0);
u.n[bottom] = neumann (0);
p[bottom] = dirichlet (0);

vector nv[];
void normal_vector (scalar s) {
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
  
  j = 0;
  Re = 1.;
  run();

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
  
}


scalar test[];
event moving_cylinder (i++) {
 
  /*
  solid (test, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
  int k;
  double vof;
  normal_vector(airfoil);
  foreach() {
    if (airfoil[] > 0 && airfoil[] < 1) {
      double lambda = fabs(nv.x[]) + fabs(nv.y[]);
      double eta = 0.065*(1 - sq(lambda)) + 0.39;
      k = 1;
      vof = airfoil[];
      // fprintf (stderr, "alpha k=%d i=%d t=%g n.x=%g n.y=%g lam=%g eta=%g vof=%g x=%g y=%g delta=%g\n",k,i,t,
      //       nv.x[], nv.y[], lambda, eta, vof, x, y, Delta);
      airfoil[] = 0.5*(1 - tanh ((sqrt(sq(x - ci.x) + sq(y - ci.y)) - (D/2))/
      				 (lambda * eta * sqrt(2) * Delta)));
      k++;
      vof = airfoil[];
      // fprintf (stderr, "alpha k=%d i=%d t=%g n.x=%g n.y=%g lam=%g eta=%g vof=%g x=%g y=%g delta=%g\n",k,i,t,
      //       nv.x[], nv.y[], lambda, eta, vof, x, y, Delta);
    }
  }
  */
}

event init (t = 0) {
  mask(y > 6 ? top: y < -6 ? bottom : none);

  solid (airfoil, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
  
  foreach()
    u.x[] = airfoil[]? 0: U0;
}


scalar airfoilsf[], airfoilsf1[], airfoilsf2[];

event properties (i++) {
  solid (airfoil, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
   
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);

  
  foreach()
    airfoilsf[] = (4.*airfoil[] + 
		   2.*(airfoil[0,1] + airfoil[0,-1] + airfoil[1,0] + airfoil[-1,0]) +
		   airfoil[-1,-1] + airfoil[1,-1] + airfoil[1,1] + airfoil[-1,1])/16.;


  foreach()
    airfoilsf1[] = (4.*airfoilsf[] + 
		    2.*(airfoilsf[0,1] + airfoilsf[0,-1] + airfoilsf[1,0] + airfoilsf[-1,0]) +
		    airfoilsf[-1,-1] + airfoilsf[1,-1] + airfoilsf[1,1] + airfoilsf[-1,1])/16.;


  foreach()
    airfoilsf2[] = (4.*airfoilsf1[] + 
		    2.*(airfoilsf1[0,1] + airfoilsf1[0,-1] + airfoilsf1[1,0] + airfoilsf1[-1,0]) +
		    airfoilsf1[-1,-1] + airfoilsf1[1,-1] + airfoilsf1[1,1] + airfoilsf1[-1,1])/16.;

  foreach()
    airfoil[] = (4.*airfoilsf2[] + 
		 2.*(airfoilsf2[0,1] + airfoilsf2[0,-1] + airfoilsf2[1,0] + airfoilsf2[-1,0]) +
		 airfoilsf2[-1,-1] + airfoilsf2[1,-1] + airfoilsf2[1,1] + airfoilsf2[-1,1])/16.;
  
  // airfoil.prolongation = refine_bilinear;
  // airfoil.dirty = true; // boundary conditions need to be updated
}

extern face vector aF;
event logfile (i++){
  coord Fp, Fmu;
  interface_force (airfoil, p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));
  
  // coord Fi = {0}, n, p;
  /*
  foreach() {
    foreach_dimension()
      n.x = nv.x[];
    coord n = interface_normal (point, airfoil);
    double alpha_p = plane_alpha (airfoil[], n);
    double area = pow(Delta, dimension - 1) * plane_area_center (n, alpha_p, &p);
    foreach_dimension()
	Fi.x += (aF.x[]*u.x[]*n.x*area);
    }
  */
  // double CD1 = (Fi.x)/(0.5*sq(U0)*(D));
  // double CL1 = (Fi.y)/(0.5*sq(U0)*(D));
  
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
  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, E);
}


event movie (t += 0.01; t <= t_end)
{
  /*
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

  char name4[80];
  sprintf (name4, "%d-stream-lines-%g", j, t);
  FILE * fp4 = fopen (name4, "w");
  view (fov = 2, tx = -0.375, ty = -0.20,
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("airfoil", "sf", filled = 1, lw = 5);
  save (fp = fp4);
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
