/*

  closepair.c

  $Revision: 1.34 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  Assumes point pattern is sorted in increasing order of x coordinate

  paircount()    count the total number of pairs (i, j) with distance < rmax

  Cclosepaircounts
                count for each i the number of j with distance < rmax

  crosscount()   count number of close pairs in two patterns

  (note: Ccrosspaircounts is defined in Estrauss.c)

  duplicatedxy() find duplicated (x,y) pairs

  Fclosepairs()  extract close pairs of coordinates 
                 .C interface - output vectors have Fixed length 

  Fcrosspairs()  extract close pairs in two patterns 
                 .C interface - output vectors have Fixed length 

  Vclosepairs()  extract close pairs of coordinates 
                 .Call interface - output vectors have Variable length 

  Vcrosspairs()  extract close pairs in two patterns 
                 .Call interface - output vectors have Variable length 

*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#define OK 0
#define ERR_OVERFLOW 1
#define ERR_ALLOC 2

#define FAILED(X) ((void *)(X) == (void *)NULL)

#define intRealloc(PTR, OLDLENGTH, NEWLENGTH) \
  (int *) S_realloc((char *) PTR, NEWLENGTH, OLDLENGTH, sizeof(int))

#define dblRealloc(PTR, OLDLENGTH, NEWLENGTH) \
  (double *) S_realloc((char *) PTR, NEWLENGTH, OLDLENGTH, sizeof(double))

double sqrt();

/* count TOTAL number of close pairs */

void paircount(nxy, x, y, rmaxi, count) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance */
     /* output */
     int *count;
{
  int n, maxchunk, i, j, counted;
  double xi, yi, rmax, r2max, dx, dy, a;

  n = *nxy;
  rmax = *rmaxi;
  r2max = rmax * rmax;

  *count = counted = 0;
  if(n == 0) 
    return;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 

  while(i < n) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n) maxchunk = n;

    for(; i < maxchunk; i++) {

      xi = x[i];
      yi = y[i];

      if(i > 0) { 
	/* scan backwards from i */
	for(j = i - 1; j >= 0; j--) {
	  dx = x[j] - xi;
	  a = r2max - dx * dx;
	  if(a < 0) 
	    break;
	  dy = y[j] - yi;
	  a -= dy * dy;
	  if(a >= 0)
	    ++counted;
	}
      }
      if(i + 1 < n) {
	/* scan forwards from i */
	for(j = i + 1; j < n; j++) {
	  dx = x[j] - xi;
	  a = r2max - dx * dx;
	  if(a < 0) 
	    break;
	  dy = y[j] - yi;
	  a -= dy * dy;
	  if(a >= 0)
	    ++counted;
	}
      } 
      /* end loop over i */
    }
  } 

  *count = counted;
}

/* count for each i the number of j closer than distance r */

void Cclosepaircounts(nxy, x, y, rmaxi, counts) 
     /* inputs */
     int *nxy;         /* number of (x,y) points */
     double *x, *y;    /* (x,y) coordinates */
     double *rmaxi;    /* maximum distance */
     /* output VECTOR, assumed initialised to 0 */
     int *counts;
{
  int n, maxchunk, i, j;
  double xi, yi, rmax, r2max, dx, dy, a;

  n = *nxy;
  rmax = *rmaxi;
  r2max = rmax * rmax;

  if(n == 0) 
    return;

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 

  while(i < n) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n) maxchunk = n;

    for(; i < maxchunk; i++) {

      xi = x[i];
      yi = y[i];

      if(i > 0) { 
	/* scan backwards from i */
	for(j = i - 1; j >= 0; j--) {
	  dx = x[j] - xi;
	  a = r2max - dx * dx;
	  if(a < 0) 
	    break;
	  dy = y[j] - yi;
	  a -= dy * dy;
	  if(a >= 0)
	    (counts[i])++;
	}
      }
      if(i + 1 < n) {
	/* scan forwards from i */
	for(j = i + 1; j < n; j++) {
	  dx = x[j] - xi;
	  a = r2max - dx * dx;
	  if(a < 0) 
	    break;
	  dy = y[j] - yi;
	  a -= dy * dy;
	  if(a >= 0)
	    (counts[i])++;
	}
      } 
      /* end loop over i */
    }
  } 
}


/*
  analogue for two different point patterns
*/


void crosscount(nn1, x1, y1, nn2, x2, y2, rmaxi, count) 
     /* inputs */
     int *nn1, *nn2;
     double *x1, *y1, *x2, *y2, *rmaxi;
     /* output */
     int *count;
{
  int n1, n2, maxchunk, i, j, jleft, counted;
  double x1i, y1i, rmax, r2max, xleft, dx, dy, a;

  n1 = *nn1;
  n2 = *nn2;
  rmax = *rmaxi;
  r2max = rmax * rmax;

  *count = counted = 0;

  if(n1 == 0 || n2 == 0) 
    return;

  jleft = 0;

  i = 0; maxchunk = 0; 

  while(i < n1) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n1) maxchunk = n1;

    for(; i < maxchunk; i++) {
  
      x1i = x1[i];
      y1i = y1[i];

      /* 
	 adjust starting index
      */
      xleft = x1i - rmax;
      while((x2[jleft] < xleft) && (jleft+1 < n2))
	++jleft;

      /* 
	 process from j=jleft until dx > rmax
      */
      for(j=jleft; j < n2; j++) {
	dx = x2[j] - x1i;
	a  = r2max - dx * dx;
	if(a < 0)
	  break;
	dy = y2[j] - y1i;
	a -= dy * dy;
	if(a > 0) 
	  ++counted;
      }
    }
  }
  *count = counted;
}



/*
  Find duplicated locations

   xx, yy are not sorted
*/


void duplicatedxy(n, x, y, out) 
     /* inputs */
     int *n;
     double *x, *y;
     /* output */
     int *out;  /* logical vector */
{
  int m, i, j;
  double xi, yi;
  m = *n;
  for(i = 1; i < m; i++) {
    R_CheckUserInterrupt();
    xi = x[i];
    yi = y[i];
    for(j = 0; j < i; j++) 
      if((x[j] == xi) && (y[j] == yi)) 
	break;
    if(j == i) out[i] = 0; else out[i] = 1;
  }
}

/* ............... fixed output length .............. */

void Fclosepairs(nxy, x, y, r, noutmax, 
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
  int n, k, kmax, maxchunk, i, j;
  double xi, yi, rmax, r2max, dx, dy, dx2, d2;

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

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < n) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n) maxchunk = n;

    for(; i < maxchunk; i++) {

      xi = x[i];
      yi = y[i];

      if(i > 0) {
	/* scan backwards */
	for(j = i - 1; j >= 0; j--) {
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 > r2max)
	    break;
	  dy = y[j] - yi;
	  d2 = dx2 + dy * dy;
	  if(d2 <= r2max) {
	    /* add this (i, j) pair to output */
	    if(k >= kmax) {
	      *nout = k;
	      *status = ERR_OVERFLOW;
	      return;
	    }
	    jout[k] = j + 1;  /* R indexing */
	    iout[k] = i + 1;
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
    
      if(i + 1 < n) {
	/* scan forwards */
	for(j = i + 1; j < n; j++) {
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 > r2max)
	    break;
	  dy = y[j] - yi;
	  d2 = dx2 + dy * dy;
	  if(d2 <= r2max) {
	    /* add this (i, j) pair to output */
	    if(k >= kmax) {
	      *nout = k;
	      *status = ERR_OVERFLOW;
	      return;
	    }
	    jout[k] = j + 1;  /* R indexing */
	    iout[k] = i + 1; 
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
  }
  *nout = k;
}

void Fcrosspairs(nn1, x1, y1, nn2, x2, y2, rmaxi, noutmax, 
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
  int n1, n2, maxchunk, k, kmax, i, j, jleft;
  double x1i, y1i, rmax, r2max, xleft, dx, dy, dx2, d2;

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

  jleft = 0;

  i = 0; maxchunk = 0; 

  while(i < n1) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > n1) maxchunk = n1;

    for(; i < maxchunk; i++) {

      x1i = x1[i];
      y1i = y1[i];

      /* 
	 adjust starting position jleft

      */
      xleft = x1i - rmax;
      while((x2[jleft] < xleft) && (jleft+1 < n2))
	++jleft;


      /* 
	 process from j=jleft until dx > rmax
      */
      for(j=jleft; j < n2; j++) {
	dx = x2[j] - x1i;
	dx2 = dx * dx;
	if(dx2 > r2max)
	  break;
	dy = y2[j] - y1i;
	d2 = dx2 + dy * dy;
	if(d2 <= r2max) {
	  /* add this (i, j) pair to output */
	  if(k >= kmax) {
	    *nout = k;
	    *status = ERR_OVERFLOW;
	    return;
	  }
	  jout[k] = j + 1;  /* R indexing */
	  iout[k] = i + 1;
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
  }
  *nout = k;
}


/* ........  versions that return variable-length vectors ......... */

#define SINGLE

/* return i, j only */
#define CLOSEFUN VcloseIJpairs
#define CROSSFUN VcrossIJpairs
#undef THRESH
#undef COORDS
#undef DIST
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS
#undef DIST

/* return i, j, d */
#define CLOSEFUN VcloseIJDpairs
#define CROSSFUN VcrossIJDpairs
#undef THRESH
#undef COORDS
#define DIST
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS
#undef DIST

/* return i, j, xi, yi, xj, yj, dx, dy, d */
#define CLOSEFUN Vclosepairs
#define CROSSFUN Vcrosspairs
#undef THRESH
#define COORDS
#define DIST
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS
#undef DIST

/* return i, j, t where t = 1{d < s} */

#define CLOSEFUN Vclosethresh
#define CROSSFUN Vcrossthresh
#define THRESH
#undef COORDS
#include "closefuns.h"
#undef CLOSEFUN
#undef CROSSFUN
#undef THRESH
#undef COORDS

