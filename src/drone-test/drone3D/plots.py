# -*- coding: utf-8 -*-
# #include "grid/octree.h"
# #include "embed.h"
# #include "navier-stokes/centered.h"
# #include "navier-stokes/perfs.h"
# #include "utils.h"
# #include "distance.h"
# #include "view.h"
# 
# int LEVEL = 8;
# double D = 4.6; // prop diameter
# double h = 0.175; // prop height
# double Re = 600000;
# double U0 = 1.;
# double omega = 10. [0, -1];// prop rotational speed
# face vector muv[];
# 
# 
# void fraction_from_stl (scalar s, face vector fv, FILE * fp) {
#   coord * p = input_stl (fp); // import CAD
#   coord min, max;
#   bounding_box (p, &min, &max);
#   double maxl = -HUGE;
#   foreach_dimension()
#     if (max.x - min.x > maxl)
#       maxl = max.x - min.x;
# 
#   // STL points to domain points
#   scalar d[];
#   distance (d, p);
#   while (adapt_wavelet ({d}, (double[]){5e-4}, LEVEL, 5).nf);
# 
#   vertex scalar phi[];
#   foreach_vertex(){
#     phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
# 	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
#   }
#   fractions (phi, s, fv);
#   fractions_cleanup (s, fv);
#   
# }
# 
# 
# int main() {
#   init_grid (32);
#   double L0 = 5;
#   size (L0);
#   origin (-L0/2, -L0/2,-L0/2);
#   mu = muv;
#   run();
# }
# 
# 
# event init (t = 0) {
# 
#   if (!restore (file = "restart")) {
#     FILE * fp = fopen ("/home/spencer/basilisk/CTFL/files/prop_test5.stl", "r");
#     if (fp == NULL)
#       fprintf(stderr, "STL file is NULL\n");
#     else {
#       fprintf(stderr, "STL is NOT null\n");
#       fraction_from_stl (cs, fs, fp);
#       foreach()
#         cs[] = 1 - cs[];
#       foreach_face()
# 	fs.x[] = 1 - fs.x[];
#       fclose (fp);
#     }
#   }
#   /*
#   FILE * cs_cad = fopen ("cs-points.txt", "w");
#   FILE * fsx_cad = fopen ("fsx-points.txt", "w");
#   FILE * fsy_cad = fopen ("fsy-points.txt", "w");
# 
#   foreach()
#     if(cs[] > 0 && cs[] < 1) {
#       fprintf (cs_cad,"%g %g %g %g\n", x, y, z, cs[]);
#       fprintf (fsx_cad,"%g %g %g %g\n", x, y, z, fs.x[]);
#       fprintf (fsy_cad,"%g %g %g %g\n", x, y, z, fs.y[]);
#     }
# 
#   fclose (cs_cad);
#   fclose (fsx_cad);
#   fclose (fsy_cad);
#   */
#   foreach() {
#     double theta = atan(y/x);
#     double r = sqrt(sq(x) + sq(y));
#     u.x[] = x > 0? r*omega*sin(theta):-r*omega*sin(theta);
#     u.y[] = x > 0? -r*omega*cos(theta): r*omega*cos(theta);
#     u.z[] = 0.;
#    }
# 
# 
# u.n[right]    = dirichlet( omega*y);
# u.n[left]     = dirichlet( omega*y);
# u.n[top]      = dirichlet(-omega*x);
# u.n[bottom]   = dirichlet(-omega*x);
# 
# 
# p[right]      = neumann (0);
# p[left]       = neumann (0);
# p[top]        = neumann (0);
# p[bottom]     = neumann (0);
# 
# 
# pf[right]     = neumann (0);
# pf[left]      = neumann (0);
# pf[top]       = neumann (0);
# pf[bottom]    = neumann (0);
# 
# 
# u.t[right]    = dirichlet(-omega*fabs(y));
# u.t[left]     = dirichlet( omega*fabs(y));
# u.t[top]      = dirichlet( omega*fabs(x));
# u.t[bottom]   = dirichlet(-omega*fabs(x));
# 
# 
# uf.n[right]   = dirichlet( omega*y);
# uf.n[left]    = dirichlet( omega*y);
# uf.n[top]     = dirichlet(-omega*x);
# uf.n[bottom]  = dirichlet(-omega*x);
# 
# uf.t[right]   = dirichlet(-omega*fabs(y));
# uf.t[left]    = dirichlet( omega*fabs(y));
# uf.t[top]     = dirichlet( omega*fabs(x));
# uf.t[bottom]  = dirichlet(-omega*fabs(x));
# 
# 
# // 3D
# u.n[front]    = dirichlet(1);
# u.t[front]    = dirichlet(1);
# p[front]      = neumann(0);
# pf[front]     = neumann(0);
# 
# u.n[back]     = neumann(0);
# u.t[back]     = neumann(0);
# p[back]       = dirichlet (0);
# pf[back]      = dirichlet (0);
# 
# u.n[embed]    = dirichlet(0);
# u.t[embed]    = dirichlet(0);
# 
# }
# 
# event properties (i++) {
#   foreach_face()
#     muv.x[] = fm.x[]*(uf.x[])*(D)/(Re);
#   boundary ((scalar *) {muv});
# }
# 
# 
# face vector av[];
# 
# event acceleration (i++) {
#   coord cor = {-2.0*omega, 2.0*omega};
#   foreach_face(x)
#     av.x[] = (cor.x*uf.y[]) + sq(omega)*x;
# 
#   foreach_face(y)
#     av.y[] = cor.y*uf.x[] + sq(omega)*y;
# 
#   foreach_face(z)
#     av.z[] = 0.;
# 
#   a = av;
# }
# 
# 
# event logfile (i++, t <= 5) {
# 
# /*
#   coord Fp, Fmu;
#   embed_force (p, u, mu, &Fp, &Fmu);
#   double CD = (Fp.x + Fmu.x)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
#   double CL = (Fp.y + Fmu.y)/(0.5*sq(omega*D/2)*3.1415*sq(D/2));
# */
#   fprintf (stderr, "%d %g\n", i, t);
# }
# 
# event movie (t += 0.01){
# }
# 
# event adapt (i++) {
#   double uemax = 0.001;
#   adapt_wavelet ({cs,u}, (double[]) {1.e-2,(0.5),(0.5),(0.5)},
# 		  maxlevel = LEVEL, minlevel = (1));
# }
