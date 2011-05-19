#include <Rmath.h>

/*

  denspt.c

  $Revision: 1.7 $     $Date: 2011/05/17 12:18:53 $

  Assumes point pattern is sorted in increasing order of x coordinate

  *denspt*     Density estimate at points
  *smoopt*     Smoothed mark values at points

*/

#define TWOPI M_2PI

double sqrt(), exp();

/* ----------------- density estimation -------------------- */

void denspt(nxy, x, y, rmaxi, sig, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *sig;      /* Gaussian sd */
     /* output */
     double *result;   /* vector of computed density values */
{
  int n, i, j, jleft, jright;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, resulti;
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

    resulti = 0.0;

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
	  resulti += exp(-d2/twosig2);
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
	  resulti += exp(-d2/twosig2);
	}
      }
    }
    /* commit */
    result[i] = coef * resulti;
  }
}


void wtdenspt(nxy, x, y, rmaxi, sig, weight, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance */
     double *sig;      /* Gaussian sd */
     double *weight;      /* vector of weights */
     /* output */
     double *result;    /* vector of weighted density values */
{
  int n, i, j, jleft, jright;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, resulti;
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

    resulti = 0.0;

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
	  resulti += weight[j] * exp(-d2/twosig2);
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
	  resulti += weight[j] * exp(-d2/twosig2);
	}
      }
    }
    result[i] = coef * resulti;
  }
}

/* ------------- anisotropic versions -------------------- */

void adenspt(nxy, x, y, rmaxi, detsigma, sinv, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result;   /* vector of density values */
{
  int n, i, j, jleft, jright;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, coef, resulti;
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

    resulti = 0.0;

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
	  resulti += exp(-(dx * (dx * s11 + dy * s12)
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
	  resulti += exp(-(dx * (dx * s11 + dy * s12)
			  + dy * (dx * s21 + dy * s22))/2.0);
	}
      }
    }
    result[i] = coef * resulti;
  }
  
}


void awtdenspt(nxy, x, y, rmaxi, detsigma, sinv, weight, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     double *weight;      /* vector of weights */
     /* output */
     double *result;    /* vector of weighted density values */
{
  int n, i, j, jleft, jright;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2, coef, resulti;
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
    resulti = 0.0;
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
	  resulti += weight[j] * exp(-(dx * (dx * s11 + dy * s12)
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
	  resulti += weight[j] * exp(-(dx * (dx * s11 + dy * s12)
					+ dy * (dx * s21 + dy * s22))/2.0);
	}
      }
    }
    result[i] = coef * resulti;
  }
}


/* --------------- smoothing --------------------------- */

void smoopt(nxy, x, y, v, self, rmaxi, sig, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *v;        /* vector of mark values to be smoothed */
     int *self;       /* 0 if leave-one-out */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *sig;      /* Gaussian sd */
     /* output */
     double *result;   /* vector of computed smoothed values */
{
  int n, i, j, jleft, jright, countself;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2;
  double sigma, twosig2;
  double numer, denom, wij; 

  n = *nxy;
  rmax = *rmaxi;
  sigma = *sig;
  countself = *self;

  twosig2 = 2.0 * sigma * sigma;

  r2max = rmax * rmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {

    numer = denom = 0.0;

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
	  wij = exp(-d2/twosig2);
	  denom += wij; 
	  numer += wij * v[j];
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
	  wij = exp(-d2/twosig2);
	  denom += wij; 
	  numer += wij * v[j];
	}
      }
    }
    /* commit */
    if(countself != 0) {
      numer += 1;
      denom += v[i];
    }
    result[i] = numer/denom;
  }
}


void wtsmoopt(nxy, x, y, v, self, rmaxi, sig, weight, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *v;        /* vector of mark values to be smoothed */
     int *self;       /* 0 if leave-one-out */
     double *rmaxi;    /* maximum distance */
     double *sig;      /* Gaussian sd */
     double *weight;      /* vector of weights */
     /* output */
     double *result;    /* vector of computed smoothed values */
{
  int n, i, j, jleft, jright, countself;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2;
  double sigma, twosig2;
  double numer, denom, wij; 

  n = *nxy;
  rmax = *rmaxi;
  sigma = *sig;
  countself = *self;

  twosig2 = 2.0 * sigma * sigma;

  r2max = rmax * rmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {

    numer = denom = 0.0;

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
	  wij = weight[j] * exp(-d2/twosig2);
	  denom += wij;
	  numer += wij * v[j];
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
	  wij = weight[j] * exp(-d2/twosig2);
	  denom += wij;
	  numer += wij * v[j];
	}
      }
    }
    if(countself != 0) {
      numer += weight[i];
      denom += weight[i] * v[i];
    }
    result[i] = numer/denom;
  }
}

/* ------------- anisotropic versions -------------------- */

void asmoopt(nxy, x, y, v, self, rmaxi, detsigma, sinv, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *v;        /* vector of mark values to be smoothed */
     int *self;       /* 0 if leave-one-out */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     /* output */
     double *result;   /* vector of smoothed values */
{
  int n, i, j, jleft, jright, countself;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2;
  double detsig, s11, s12, s21, s22;
  double numer, denom, wij; 

  n = *nxy;
  rmax = *rmaxi;
  detsig = *detsigma;
  countself = *self;

  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  r2max = rmax * rmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {

    numer = denom = 0.0;

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
	  wij = exp(-(dx * (dx * s11 + dy * s12)
		      + dy * (dx * s21 + dy * s22))/2.0);
	  denom += wij;
	  numer += wij * v[j];
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
	  wij = exp(-(dx * (dx * s11 + dy * s12)
		      + dy * (dx * s21 + dy * s22))/2.0);
	  denom += wij;
	  numer += wij * v[j];
	}
      }
    }
    if(countself != 0) {
      numer += 1;
      denom += v[i];
    }
    result[i] = numer/denom;
  }
  
}


void awtsmoopt(nxy, x, y, v, self, rmaxi, detsigma, sinv, weight, result) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *v;        /* vector of mark values to be smoothed */
     int *self;       /* 0 if leave-one-out */
     double *rmaxi;    /* maximum distance at which points contribute */
     double *detsigma;  /* determinant of variance matrix */
     double *sinv;      /* inverse variance matrix (2x2, flattened) */
     double *weight;      /* vector of weights */
     /* output */
     double *result;    /* vector of smoothed values */
{
  int n, i, j, jleft, jright, countself;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2;
  double detsig, s11, s12, s21, s22;
  double numer, denom, wij; 

  n = *nxy;
  rmax = *rmaxi;
  detsig = *detsigma;
  countself = *self;

  s11 = sinv[0];
  s12 = sinv[1];
  s21 = sinv[2];
  s22 = sinv[3];

  r2max = rmax * rmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {
    numer = denom = 0.0;
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
	  wij = weight[j] * exp(-(dx * (dx * s11 + dy * s12)
				  + dy * (dx * s21 + dy * s22))/2.0);
	  denom += wij;
	  numer += wij * v[j];
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
	  wij = weight[j] * exp(-(dx * (dx * s11 + dy * s12)
				  + dy * (dx * s21 + dy * s22))/2.0);
	  denom += wij;
	  numer += wij * v[j];
	}
      }
    }
    if(countself != 0) {
      numer += weight[i];
      denom += weight[i] * v[i];
    }
    result[i] = numer/denom;
  }
}


