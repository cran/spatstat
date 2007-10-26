/*

  closepair.c

  $Revision: 1.10 $     $Date: 2007/10/26 14:55:32 $

  Assumes point pattern is sorted in increasing order of x coordinate


  paircount()    simply count the number of pairs (i, j) with distance < rmax

  closepairs()   extract all pairs of coordinates with distance < rmax

*/

#define OK 0
#define OVERFLOW 1

double sqrt();

void paircount(nxy, x, y, rmaxi, count) 
     int *nxy;
     double *x, *y, *rmaxi;
     int *count;
{
  int n, i, j, jleft, jright, counted;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2;

  n = *nxy;
  rmax = *rmaxi;
  r2max = rmax * rmax;

  *count = counted = 0;
  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {
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
	if(d2 <= r2max) 
	  ++counted;
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
	if(d2 <= r2max) 
	  ++counted;
      }
    }
  }
  *count = counted;
}


void closepairs(nxy, x, y, r, noutmax, 
	      nout, iout, jout, 
	      xiout, yiout, xjout, yjout, dxout, dyout, dout,
	      status)
     /* inputs */
     int *nxy, *noutmax;
     double *x, *y, *r;
     /* outputs */
     int *nout, *iout, *jout;
     double *xiout, *yiout, *xjout, *yjout, *dxout, *dyout, *dout;
     int *status;
{
  int n, k, kmax, i, j, jleft, jright;
  double xi, yi, rmax, r2max, xleft, xright, dx, dy, d2;

  n = *nxy;
  rmax = *r;
  r2max = rmax * rmax;

  *status = OK;
  *nout = 0;
  k = 0;   /* k is the next available storage location 
              and also the current length of the list */ 
  kmax = *noutmax;

  if(n == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n; i++) {
    xi = x[i];
    yi = y[i];
    /* search all points with x in [xleft, xright] */
    xleft = xi - rmax;
    xright = xi + rmax;

    /* 
       adjust scope of search [jleft, jright]

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
	  /* add this (i, j) pair to output */
	  if(k >= kmax) {
	    *nout = k;
	    *status = OVERFLOW;
	    return;
	  }
	  jout[k] = j;
	  iout[k] = i;
	  xiout[k] = xi;
	  yiout[k] = yi;
	  xjout[k] = x[j];
	  yjout[k] = y[j];
	  dxout[k] = dx;
	  dyout[k] = dy;
	  dout[k] = sqrt(d2);
	  ++k;
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
	if(d2 <= r2max) {
	  /* add this (i, j) pair to output */
	  if(k >= kmax) {
	    *nout = k;
	    *status = OVERFLOW;
	    return;
	  }
	  jout[k] = j;
	  iout[k] = i;
	  xiout[k] = xi;
	  yiout[k] = yi;
	  xjout[k] = x[j];
	  yjout[k] = y[j];
	  dxout[k] = dx;
	  dyout[k] = dy;
	  dout[k] = sqrt(d2);
	  ++k;
	}
      }
    }
  }
  *nout = k;
}

/*
  analogue for two different point patterns
*/

void crosscount(nn1, x1, y1, nn2, x2, y2, rmaxi, count) 
     int *nn1, *nn2;
     double *x1, *y1, *x2, *y2, *rmaxi;
     int *count;
{
  int n1, n2, i, j, jleft, jright, counted;
  double x1i, y1i, rmax, r2max, xleft, xright, dx, dy, d2;

  n1 = *nn1;
  n2 = *nn2;
  rmax = *rmaxi;
  r2max = rmax * rmax;

  *count = counted = 0;

  if(n1 == 0 || n2 == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n1; i++) {
    x1i = x1[i];
    y1i = y1[i];
    /* search all points with x in [xleft, xright] */
    xleft = x1i - rmax;
    xright = x1i + rmax;

    /* 
       adjust scope of search [jleft, jright]

    */
    while((jleft+1 < n2) && x2[jleft] < xleft)
      ++jleft;

    while((jright+1 < n2) && x2[jright+1] <= xright)
      ++jright;

    /* 
       process from jleft to jright
    */
    for(j=jleft; j <= jright; j++) {
      /* squared interpoint distance */
      dx = x2[j] - x1i;
      dy = y2[j] - y1i;
      d2= dx * dx + dy * dy;
      if(d2 <= r2max) 
	++counted;
    }
  }
  *count = counted;
}


void crosspairs(nn1, x1, y1, nn2, x2, y2, rmaxi, noutmax, 
	      nout, iout, jout, 
	      xiout, yiout, xjout, yjout, dxout, dyout, dout,
	      status)
     /* inputs */
     int *nn1, *nn2, *noutmax;
     double *x1, *y1, *x2, *y2, *rmaxi;
     /* outputs */
     int *nout, *iout, *jout;
     double *xiout, *yiout, *xjout, *yjout, *dxout, *dyout, *dout;
     int *status;
{
  int n1, n2, k, kmax, i, j, jleft, jright;
  double x1i, y1i, rmax, r2max, xleft, xright, dx, dy, d2;

  n1 = *nn1;
  n2 = *nn2;
  rmax = *rmaxi;
  r2max = rmax * rmax;

  *status = OK;
  *nout = 0;
  k = 0;   /* k is the next available storage location 
              and also the current length of the list */ 
  kmax = *noutmax;

  if(n1 == 0 || n2 == 0) 
    return;

  jleft = jright = 0;

  for(i = 0; i < n1; i++) {
    x1i = x1[i];
    y1i = y1[i];
    /* search all points with x in [xleft, xright] */
    xleft = x1i - rmax;
    xright = x1i + rmax;

    /* 
       adjust scope of search [jleft, jright]

    */
    while((jleft+1 < n2) && x2[jleft] < xleft)
      ++jleft;

    while((jright+1 < n2) && x2[jright+1] <= xright)
      ++jright;

    /* 
       process from jleft to jright
    */
    for(j=jleft; j <= jright; j++) {
      /* squared interpoint distance */
      dx = x2[j] - x1i;
      dy = y2[j] - y1i;
      d2= dx * dx + dy * dy;
      if(d2 <= r2max) {
	/* add this (i, j) pair to output */
	if(k >= kmax) {
	  *nout = k;
	  *status = OVERFLOW;
	  return;
	}
	jout[k] = j;
	iout[k] = i;
	xiout[k] = x1i;
	yiout[k] = y1i;
	xjout[k] = x2[j];
	yjout[k] = y2[j];
	dxout[k] = dx;
	dyout[k] = dy;
	dout[k] = sqrt(d2);
	++k;
      }
    }
  }
  *nout = k;
}

/*
  Find duplicated locations

   xx, yy are not sorted
*/


void duplicatedxy(n, x, y, out) 
     int *n;
     double *x, *y;
     int *out;
{
  int m, i, j;
  double xi, yi;
  m = *n;
  for(i = 1; i < m; i++) {
    xi = x[i];
    yi = y[i];
    for(j = 0; j < i; j++) 
      if((x[j] == xi) && (y[j] == yi)) 
	break;
    if(j == i) out[i] = 0; else out[i] = 1;
  }
}


