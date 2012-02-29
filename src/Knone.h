/*
  
  Knone.h

  Code template for K function estimators in Knone.c

  Variables:

     FNAME        function name

     OUTTYPE      storage type of the output 'numer' 
                  ('int' or 'double')

     WEIGHTED     #defined for weighted (inhom) K function


  $Revision: 1.2 $     $Date: 2012/02/29 04:57:21 $

*/

void FNAME(
	   nxy, x, y, 
#ifdef WEIGHTED
	   w,
#endif
	   nr, rmax, numer) 
     /* inputs */
     int *nxy, *nr;
     double *x, *y, *rmax;
#ifdef WEIGHTED
     double *w;
#endif
     /* output */
     OUTTYPE *numer;
{
  int i, j, l, n, nt, n1, lmin, lmax;
  double dt, tmax, tmax2, xi, yi;
  double xleft, xright, dratio, dij, dij2, dx, dy;
#ifdef WEIGHTED
  double wi, wj, wij;
#endif

#ifdef WEIGHTED

#define ZERO 0.0
#define WI wi
#define WJ wj
#define WIJ wij

#else 

#define ZERO 0
#define WI 1
#define WJ 1
#define WIJ 1

#endif

  n = *nxy;
  nt = *nr;

  n1 = n - 1;
  lmax = nt - 1;

  dt = (*rmax)/(nt-1);
  tmax = *rmax;
  tmax2 = tmax * tmax;

  /* initialise */
  for(l = 0; l < nt; l++)
    numer[l] = ZERO;

  if(n == 0) 
    return;

  for(i = 0; i < n; i++) {
#ifdef WEIGHTED
    wi = w[i];
#endif
    xi = x[i];
    yi = y[i];

    /* scan through points (x[j],y[j]) */
    xleft = xi - tmax;
    xright = xi + tmax;

    /* 
       scan backward from i-1 
       until x < xleft or until we run out 
    */
    if(i > 0) {
      for(j=i-1; x[j] >= xleft && j >= 0; j--) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	dij2= dx * dx + dy * dy;
	if(dij2 < tmax2) {
#ifdef WEIGHTED 
	  wj = w[j];
#endif
	  /* increment numerator for all r >= dij */
	  dij = (double) sqrt(dij2);
	  dratio = dij/dt;
	  /* smallest integer greater than or equal to dratio */
	  lmin = (int) ceil(dratio);
	  /* increment entries lmin to lmax inclusive */
	  if(lmin <= lmax) {
#ifdef WEIGHTED
	    wij = wi * wj;
#endif
	    for(l = lmin; l <= lmax; l++) 
	      numer[l] += WIJ;
	  }
	}
      }
    }

    /* 
       scan forward from i+1 
       until x > xright or until we run out 
    */
    if(i < n1) {
      for(j=i+1; x[j] <= xright && j < n; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	dij2= dx * dx + dy * dy;
	if(dij2 < tmax2) {
#ifdef WEIGHTED 
	  wj = w[j];
#endif
	  /* increment numerator for all r >= dij */
	  dij = (double) sqrt(dij2);
	  dratio = dij/dt;
	  /* smallest integer greater than or equal to dratio */
	  lmin = (int) ceil(dratio);
	  /* increment entries lmin to lmax inclusive */
	  if(lmin <= lmax) {
#ifdef WEIGHTED
	    wij = wi * wj;
#endif
	    for(l = lmin; l <= lmax; l++) 
	      numer[l] += WIJ;
	  }
	}
      }
    }
  }
}

#undef ZERO
#undef WI 
#undef WJ 
#undef WIJ

