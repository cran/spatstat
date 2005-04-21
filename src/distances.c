/*

  distances.c

  Distances between points

  $Revision: 1.8 $     $Date: 2005/04/21 02:42:29 $

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


/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

void nndistsort(n, x, y, nnd, huge)
     int *n;
     double *x, *y, *nnd, *huge;
{ 
  int npoints, i, left, right;
  double d, dmin, d2, d2min, xi, yi, dx, dy, hu, hu2;
  /* FILE *out; */

  hu = *huge;
  hu2 = hu * hu;

  /* out = fopen("out", "w"); */
  npoints = *n;

  for(i = 0; i < npoints; i++) {
    /* fprintf(out, "\ni=%d\n", i); */
    dmin = hu;
    d2min = hu2;
    xi = x[i];
    yi = y[i];
    /* search backward */
    for(left = i - 1;
        left >= 0 && (dy = (yi - y[left])) < dmin ;
	--left)
      {
	/* fprintf(out, "L"); */
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
	/* fprintf(out, "R"); */
	dx = x[right] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	}
      }
    /* fprintf(out, "\n"); */
    nnd[i] = dmin;
  }
  /* fclose(out); */
}
