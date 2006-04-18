/*

  distances.c

  Distances between points

  $Revision: 1.12 $     $Date: 2006/04/15 17:11:45 $

 */

#include <math.h>
/* #include <stdio.h> */

double sqrt();

void pairdist(n, x, y, d)
     int *n;
     double *x, *y, *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, dx, dy, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dist = sqrt( dx * dx + dy * dy ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

/* squared distances */

void pair2dist(n, x, y, d)
     int *n;
     double *x, *y, *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, dx, dy, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dist = dx * dx + dy * dy; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

void crossdist(nfrom, xfrom, yfrom, nto, xto, yto, d)
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto, *d;
{ 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, dx, dy;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	*dptr = sqrt( dx * dx + dy * dy ); 
    }
  }
}

/* squared distances */

void cross2dist(nfrom, xfrom, yfrom, nto, xto, yto, d)
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto, *d;
{ 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, dx, dy;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	*dptr = dx * dx + dy * dy; 
    }
  }
}


/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

void nndistsort(n, x, y, nnd, huge)
     int *n;
     double *x, *y, *nnd, *huge;
{ 
  int npoints, i, left, right;
  double d, dmin, d2, d2min, xi, yi, dx, dy, hu, hu2;

#ifdef SPATSTAT_DEBUG
  FILE *out; 
#endif

  hu = *huge;
  hu2 = hu * hu;

#ifdef SPATSTAT_DEBUG
  out = fopen("out", "w");
#endif

  npoints = *n;

  for(i = 0; i < npoints; i++) {

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\ni=%d\n", i); 
#endif

    dmin = hu;
    d2min = hu2;
    xi = x[i];
    yi = y[i];
    /* search backward */
    for(left = i - 1;
        left >= 0 && (dy = (yi - y[left])) < dmin ;
	--left)
      {

#ifdef SPATSTAT_DEBUG
	fprintf(out, "L");
#endif

	dx = x[left] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dy = (y[right] - yi)) < dmin ;
	++right)
      {

#ifdef SPATSTAT_DEBUG
	fprintf(out, "R");
#endif
	dx = x[right] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	}
      }
#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n");
#endif

    nnd[i] = dmin;
  }

#ifdef SPATSTAT_DEBUG
  fclose(out);
#endif

}

/* nnwhichsort: same as nndistsort, 
   but also returns id of nearest neighbour 
*/

void nnwhichsort(n, x, y, nnd, nnwhich, huge)
     int *n, *nnwhich;
     double *x, *y, *nnd, *huge;
{ 
  int npoints, i, left, right, which;
  double d, dmin, d2, d2min, xi, yi, dx, dy, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  for(i = 0; i < npoints; i++) {
    dmin = hu;
    d2min = hu2;
    which = -1;
    xi = x[i];
    yi = y[i];
    /* search backward */
    for(left = i - 1;
        left >= 0 && (dy = (yi - y[left])) < dmin ;
	--left)
      {
	dx = x[left] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	  which = left;
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dy = (y[right] - yi)) < dmin ;
	++right)
      {
	dx = x[right] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	  which = right;
	}
      }
    nnd[i] = dmin;
    nnwhich[i] = which;
  }
}

