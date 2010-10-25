#include <Rmath.h>

/*

  denspt.c

  $Revision: 1.3 $     $Date: 2010/10/24 13:18:50 $

  Assumes point pattern is sorted in increasing order of x coordinate

  Density estimate at points

*/

#define TWOPI M_2PI

double sqrt(), exp();

void denspt(nxy, x, y, rmaxi, sig, value) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *sig;      /* Gaussian sd */
     /* output */
     double *value;   /* vector of density values */
{
  int n, i, j, jleft, jright, counted;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, valuei;
  double sigma, coef, twosig2;

  n = *nxy;
  rmax = *rmaxi;
  sigma = *sig;

  twosig2 = 2.0 * sigma * sigma;
  coef = 1.0/(TWOPI * sigma * sigma);

  r2max = rmax * rmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {

    valuei = 0.0;

    xi = x[i];
    yi = y[i];

    /* search all points with x in [xleft, xright] */
    xleft = xi - rmax;
    xright = xi + rmax;

    /* 
       adjust moving boundary - scope of search [jleft, jright]

    */
    while(x[jleft] < xleft && jleft < i)
      ++jleft;

    while(jright+1 < n && x[jright+1] <= xright)
      ++jright;

    /* 
       process from jleft to i-1 
    */
    if(jleft < i) {
      for(j=jleft; j < i; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2 <= r2max) {
	  valuei += exp(-d2/twosig2);
	}
      }
    }

    /* 
       process from i+1 to jright
    */
    if(jright > i) {
      for(j=i+1; j <= jright; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2<= r2max){
	  valuei += exp(-d2/twosig2);
	}
      }
    }
    /* commit */
    value[i] = coef * valuei;
  }
}


void wtdenspt(nxy, x, y, rmaxi, sig, weight, value) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance */
     double *sig;      /* Gaussian sd */
     double *weight;      /* vector of weights */
     /* output */
     double *value;    /* vector of weighted density values */
{
  int n, i, j, jleft, jright, counted;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, valuei;
  double sigma, coef, twosig2;

  n = *nxy;
  rmax = *rmaxi;
  sigma = *sig;

  twosig2 = 2.0 * sigma * sigma;
  coef = 1.0/(TWOPI * sigma * sigma);

  r2max = rmax * rmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {

    valuei = 0.0;

    xi = x[i];
    yi = y[i];

    /* search all points with x in [xleft, xright] */
    xleft = xi - rmax;
    xright = xi + rmax;

    /* 
       adjust moving boundary - scope of search [jleft, jright]

    */
    while(x[jleft] < xleft && jleft < i)
      ++jleft;

    while(jright+1 < n && x[jright+1] <= xright)
      ++jright;

    /* 
       process from jleft to i-1 
    */
    if(jleft < i) {
      for(j=jleft; j < i; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2 <= r2max) {
	  valuei += weight[j] * exp(-d2/twosig2);
	}
      }
    }

    /* 
       process from i+1 to jright
    */
    if(jright > i) {
      for(j=i+1; j <= jright; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2<= r2max){
	  valuei += weight[j] * exp(-d2/twosig2);
	}
      }
    }
    value[i] = coef * valuei;
  }
}

/* ------------- anisotropic versions -------------------- */

void adenspt(nxy, x, y, rmaxi, detsigma, sinv, value) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *value;   /* vector of density values */
{
  int n, i, j, jleft, jright, counted;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, coef, valuei;
  double detsig, s11, s12, s21, s22;

  n = *nxy;
  rmax = *rmaxi;
  detsig = *detsigma;
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  r2max = rmax * rmax;
  coef = 1.0/(TWOPI * sqrt(detsig));

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {

    valuei = 0.0;

    xi = x[i];
    yi = y[i];
    /* search all points with x in [xleft, xright] */
    xleft = xi - rmax;
    xright = xi + rmax;

    /* 
       adjust moving boundary - scope of search [jleft, jright]

    */
    while(x[jleft] < xleft && jleft < i)
      ++jleft;

    while(jright+1 < n && x[jright+1] <= xright)
      ++jright;

    /* 
       process from jleft to i-1 
    */
    if(jleft < i) {
      for(j=jleft; j < i; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2 <= r2max) {
	  valuei += exp(-(dx * (dx * s11 + dy * s12)
			  + dy * (dx * s21 + dy * s22))/2.0);
	}
      }
    }

    /* 
       process from i+1 to jright
    */
    if(jright > i) {
      for(j=i+1; j <= jright; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2<= r2max){
	  valuei += exp(-(dx * (dx * s11 + dy * s12)
			  + dy * (dx * s21 + dy * s22))/2.0);
	}
      }
    }
    value[i] = coef * valuei;
  }
  
}


void awtdenspt(nxy, x, y, rmaxi, detsigma, sinv, weight, value) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     double *weight;      /* vector of weights */
     /* output */
     double *value;    /* vector of weighted density values */
{
  int n, i, j, jleft, jright, counted;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, coef, valuei;
  double detsig, s11, s12, s21, s22;

  n = *nxy;
  rmax = *rmaxi;
  detsig = *detsigma;
  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  r2max = rmax * rmax;
  coef = 1.0/(TWOPI * sqrt(detsig));

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {
    valuei = 0.0;
    xi = x[i];
    yi = y[i];
    /* search all points with x in [xleft, xright] */
    xleft = xi - rmax;
    xright = xi + rmax;

    /* 
       adjust moving boundary - scope of search [jleft, jright]

    */
    while(x[jleft] < xleft && jleft < i)
      ++jleft;

    while(jright+1 < n && x[jright+1] <= xright)
      ++jright;

    /* 
       process from jleft to i-1 
    */
    if(jleft < i) {
      for(j=jleft; j < i; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2 <= r2max) {
	  valuei += weight[j] * exp(-(dx * (dx * s11 + dy * s12)
					+ dy * (dx * s21 + dy * s22))/2.0);
	}
      }
    }

    /* 
       process from i+1 to jright
    */
    if(jright > i) {
      for(j=i+1; j <= jright; j++) {
	/* squared interpoint distance */
	dx = x[j] - xi;
	dy = y[j] - yi;
	d2= dx * dx + dy * dy;
	if(d2<= r2max){
	  valuei += weight[j] * exp(-(dx * (dx * s11 + dy * s12)
					+ dy * (dx * s21 + dy * s22))/2.0);
	}
      }
    }
    value[i] = coef * valuei;
  }
}


