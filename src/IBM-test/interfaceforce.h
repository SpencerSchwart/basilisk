#ifndef BASILISK_HEADER_24
#define BASILISK_HEADER_24
#line 1 "./../interfaceforce.h"
extern scalar airfoil;
extern face vector sf;

void boundary_cells (scalar c, scalar ink) {

  scalar phi[];

  foreach() {
    if(c[] > 1e-6 && c[] < 1. - 1e-6)
      phi[] = 1.;
    else
      phi[] = -1.;
  }
  boundary ({phi});
  fractions (phi, ink);
}


#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradiento_x (Point point, scalar s, scalar cs,
					   coord n, coord p, double bc,
					   double * coef)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !sf.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (sf.x[i + (i < 0),j] && sf.y[i,j] && sf.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = sf.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!sf.y[i,j,k+m] || !sf.y[i,j+1,k+m] ||
	    !sf.z[i,j+m,k] || !sf.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }
  *coef = 0.;
  if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}


double dirichlet_gradiento (Point point, scalar s, scalar cs,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradiento_x (point, s, cs, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradiento_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradiento_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradiento_z (point, s, cs, n, p, bc, coef);
#endif // dimension == 3
  return nodata;
}


double embed_geometryo (Point point, coord * b, coord * n)
{
  *n = facet_normal (point, airfoil, sf);
  double alpha = plane_alpha (airfoil[], *n);
  double area = plane_area_center (*n, alpha, b);
  normalize (n);
  return area;
}

double embed_interpolateo (Point point, scalar s, coord p)
{
  assert (dimension == 2);
  int i = sign(p.x), j = sign(p.y);
  if (airfoil[i] && airfoil[0,j] && airfoil[i,j])
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (airfoil[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (airfoil[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}


static inline
coord embed_gradiento (Point point, vector u, coord p, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet = false;
    double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
    // double vb = 1.;
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradiento (point, u.x, airfoil, n, p, vb, &val);
    }
    else // Neumann
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}


double interface_force (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu)
{
  double Fn = 0.;
  coord Fps = {0}, Fmus = {0};
  foreach ()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n, b;
      double area = embed_geometryo (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      Fn = area*embed_interpolateo (point, p, b);
      foreach_dimension()
	Fps.x -= Fn*n.x;

      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fm.x[] + fm.x[1];
	}
	mua /= fa;
	assert (dimension == 2);
	coord dudn = embed_gradiento (point, u, b, n);
	foreach_dimension()
	  Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
      }
    }

  *Fp = Fps; *Fmu = Fmus;
}

double embed_vorticityo (Point point, vector u, coord p, coord n)
{
  coord dudn = embed_gradient (point, u, p, n);

  return dudn.y*n.x - dudn.x*n.y;
}


#endif
