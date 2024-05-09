struct sOutputSlice {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  coord plane;	
};

trace
void output_slice (struct sOutputSlice p)
{
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  if (!p.plane.x) p.plane.x = 1.;
  if (!p.plane.y) p.plane.y = 1.;
  if (!p.plane.z) p.plane.z = 0.;
p.n++;
  
  int len = list_len(p.list);
  double ** field = (double **) matrix_new (p.n, p.n, len*sizeof(double));
  double Delta = 0.999999*L0/(p.n - 1);

  for (int i = 0; i < p.n; i++) {
    double varCoord1 = Delta*i;  // some clever way of implementing general variation of coordinates instead of mapping them directly to x, y or z
    bool varX = !(p.plane.x < 1.);
    double x = (!varX ? p.plane.x*L0 : varCoord1) + X0;
    
    for (int j = 0; j < p.n; j++) {
      double varCoord2 = Delta*j; 
      double y = (varX ? (p.plane.y < 1. ? p.plane.y*L0 : varCoord2) : varCoord1) + Y0;
      double z = (p.plane.z < 1. ? p.plane.z*L0 : varCoord2) + Z0;
      if (p.linear) {
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = interpolate (s, x, y, z);
      }
      else {
	Point point = locate (x, y, z);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif

    fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
    int k = 0; 
    for (scalar s in p.list) {
    fprintf (p.fp, "%s\n", s.name);	 
    for (int i = 0; i < p.n; i++) {
      for (int j = 0; j < p.n; j++) {
        fprintf (p.fp, "%g\t", (float) field[i][len*j + k]);	
      }
      fputc ('\n', p.fp);
    }
    k++;
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif

  matrix_free (field);
}

