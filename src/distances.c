/*

  distances.c

  Distances between pairs of points

  $Revision: 1.24 $     $Date: 2010/01/05 03:30:26 $

  pairdist      Pairwise distances
  pair2dist     Pairwise distances squared
  pairPdist     Pairwise distances with periodic correction
  pairP2dist    Pairwise distances squared, with periodic correction

  crossdist     Pairwise distances for two sets of points
  cross2dist    Pairwise distances squared, for two sets of points
  crossPdist    Pairwise distances for two sets of points, periodic correction

  matchxy       Find matches between two sets of points   

 */

#include <math.h>
/* #include <stdio.h> */

double sqrt();

void pairdist(n, x, y, d)
     /* inputs */
     int *n;
     double *x, *y;
     /* output */
     double *d;
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
     /* inputs */
     int *n;
     double *x, *y;
     /* output */
     double *d;
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
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto;
     /* output */
     double *d;
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
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto;
     /* output */
     double *d;
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



/* distances with periodic correction */

void pairPdist(n, x, y, xwidth, yheight, d)
     /* inputs */
     int *n;
     double *x, *y, *xwidth, *yheight;
     /* output */
     double *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, dx, dy, dx2, dy2, dx2p, dy2p, dist, wide, high;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;

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
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dist = sqrt( dx2p + dy2p ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

/* same function without the sqrt */

void pairP2dist(n, x, y, xwidth, yheight, d)
     /* inputs */
     int *n;
     double *x, *y, *xwidth, *yheight;
     /* output */
     double *d;
{ 
  int i, j, npoints; 
  double *dp;
  double xi, yi, dx, dy, dx2, dy2, dx2p, dy2p, dist, wide, high;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;

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
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dist = dx2p + dy2p; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

void crossPdist(nfrom, xfrom, yfrom, nto, xto, yto, xwidth, yheight, d)
     /* inputs */
     int *nto, *nfrom;
     double *xfrom, *yfrom, *xto, *yto, *xwidth, *yheight;
     /* output */
     double *d;
{ 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, dx, dy, dx2, dy2, dx2p, dy2p, wide, high;

  nf = *nfrom;
  nt = *nto;
  wide = *xwidth;
  high = *yheight;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  *dptr = sqrt( dx2p + dy2p ); 
    }
  }
}

/*

  matchxy

  Find matches between two lists of points

 */

void matchxy(na, xa, ya, nb, xb, yb, match)
     /* inputs */
     int *na, *nb;
     double *xa, *ya, *xb, *yb;
     /* output */
     int *match; 
{ 
  int i, j, Na, Nb; 
  double xai, yai;

  Na = *na;
  Nb = *nb;

  for (i=1; i < Na; i++) 
    {
      xai = xa[i];
      yai = ya[i];
      match[i] = 0;
      for (j=0; j < Nb; j++) 
	if(xai == xb[j] && yai == yb[j]) {
	  match[i] = j;
	  break;
	}
    }
}

