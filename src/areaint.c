#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

#define NGRID 16

/*
  Conditional intensity function for an area-interaction process:

  cif = eta^(1-B) where B = (uncovered area)/(pi r^2)

*/

/* Format for storage of parameters and precomputed/auxiliary data */

typedef struct AreaInt {
  /* model parameters */
  double eta;
  double r;
  /* transformations of the parameters */
  double r2;
  double range2;
  double logeta;
  double hard;
  /* periodic distance */
  double *period;
  int per;
  /* grid counting */
  double dx;
  double xgrid0;
  int *my;
  int kdisc;
  /* auxiliary indicators */
  char *neigh;
} AreaInt;

/* initialiser function */

Cdata *areaintInit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  double r, dx, dy, x0;
  int i, my, kdisc;
  AreaInt *areaint;

  /* create storage */
  areaint = (AreaInt *) R_alloc(1, sizeof(AreaInt));
  /* Interpret model parameters*/
  areaint->eta    = model.ipar[0];
  areaint->r      = r = model.ipar[1]; 
  areaint->r2     = r * r;
  areaint->range2 = 4 * r * r;    /* square of interaction distance */
  /* is the model numerically equivalent to hard core ? */
  areaint->hard   = (areaint->eta < DOUBLE_EPS);
  areaint->logeta = (areaint->hard) ? 0 : log(areaint->eta);
  /* periodic boundary conditions? */
  areaint->period = model.period;
  areaint->per    = (model.period[0] > 0.0);
  /* grid counting */
  dx = dy = areaint->dx = (2 * r)/NGRID;
  areaint->xgrid0 = -r + dx/2;
  areaint->my = (int *) R_alloc((long) NGRID, sizeof(int));
  kdisc = 0;
  for(i = 0; i < NGRID; i++) {
    x0 = areaint->xgrid0 + i * dx;
    my = floor(sqrt(r * r - x0 * x0)/dy);
    my = (my < 0) ? 0 : my;
    areaint->my[i] = my;
    kdisc += 2 * my + 1;
  }
  areaint->kdisc = kdisc;
  /* allocate space for neighbour indicators */
  areaint->neigh = (char*) R_alloc((long) state.npmax, sizeof(char));
  return((Cdata *) areaint);
}

/* conditional intensity evaluator */

double areaintCif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, kount, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double r2, d2, dx, dy, range2;
  double a, xgrid0, ygrid, ygrid0, xgrid, covfrac, cifval;
  int kdisc, kx, my, ky, covered;
  AreaInt *areaint;

  areaint = (AreaInt *) cdata;

  period = areaint->period;
  r2     = areaint->r2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;
  if(npts == 0) return ((double) 1.0);

  r2 = areaint->r2;
  dy = dx = areaint->dx;
  range2 = areaint->r2;    /* square of interaction distance */

  if(!areaint->per) {
    /*
      Euclidean distance
      First identify which data points are neighbours of (u,v)
    */
    ixp1 = ix + 1;
    /* If ix = NONE = -1, then ixp1 = 0 is correct */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	a = range2 - pow(u - x[j], 2);
	areaint->neigh[j] = FALSE;
	if(a > 0.) {
	  a -= pow(v - y[j], 2);
	  if(a > 0.)
	    areaint->neigh[j] = TRUE;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j < npts; j++) {
	a = range2 - pow(u - x[j], 2);
	areaint->neigh[j] = FALSE;
	if(a > 0.) {
	  a -= pow(v - y[j], 2);
	  if(a > 0.)
	    areaint->neigh[j] = TRUE;
	}
      }
    }
    /* scan a grid of points centred at (u,v) */
    kount = 0;
    xgrid0 = u + areaint->xgrid0;
    for(kx=0; kx<NGRID; kx++) {
      xgrid = xgrid0 + kx * dx;
      my = areaint->my[kx];
      for(ky=(-my); ky<=my; ky++) {
	ygrid = v + ky * dy;
	/*
	  Grid point (xgrid,ygrid) is inside disc of
	  radius r centred at (u,v)

	  Loop through all data points to determine
	  whether the grid point is covered by another disc
	*/
	covered = FALSE;
	if(ix > 0) {
	  for(j=0; j < ix; j++) {
	    if(areaint->neigh[j]) {
	      a = r2 - pow(xgrid - x[j], 2);
	      if(a > 0) {
		a -= pow(ygrid - y[j], 2);
		if(a > 0) {
		  /* point j covers grid point */
		  covered = TRUE;
		}
	      }
	    }
	  }
	}
	if(!covered && ixp1 < npts) {
	  for(j=ixp1; j<npts; j++) {
	    if(areaint->neigh[j]) {
	      a = r2 - pow(xgrid - x[j], 2);
	      if(a > 0) {
		a -= pow(ygrid - y[j], 2);
		if(a > 0) {
		  /* point j covers grid point */
		  covered = TRUE;
		}
	      }
	    }
	  }
	}
	/* finished scanning all data points j  */
	if(!covered)
	  kount = kount + 1;
	/* finished consideration of grid point (xgrid, ygrid) */
      }
    }
    kdisc = areaint->kdisc;
  } else {
    /*
      periodic distance
      First identify which data points are neighbours of (u,v)
    */
    ixp1 = ix + 1;
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],period);
	areaint->neigh[j] = (d2 < range2);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],period);
	areaint->neigh[j] = (d2 < range2);
      }
    }
    /* scan a grid of ngrid * ngrid points centred at (u,v) */
    kount = 0;
    kdisc = 0;
    xgrid0 = u + areaint->xgrid0;
    ygrid0 = u + areaint->xgrid0;
    for(kx=0; kx<NGRID; kx++) {
      xgrid = xgrid0 + kx * dx;
      for(ky=0; ky<NGRID; ky++) {
	ygrid = ygrid0 + ky * dy;
	/* Determine whether grid point is inside disc */
	d2 = dist2(u,v,xgrid,ygrid,period);
	if(d2 < r2) {
	  /* Grid point is inside disc of radius r centred at (u,v) */
	  kdisc = kdisc + 1;
	  /*
	    Loop through all data points to determine
	    whether the grid point is covered by another disc
	  */
	  covered = FALSE;
	  if(ix > 0) {
	    for(j=0; j< ix; j++) {
	      if(areaint->neigh[j]) {
		d2 = dist2(xgrid,ygrid,x[j],y[j],period);
		if(d2 < r2) {
		  /* point j covers grid point */
		  covered = TRUE;
		}
	      }
	    }
	  }
	  if(!covered && ixp1 < npts) {
	    for(j=ixp1; j<npts; j++) {
	      if(areaint->neigh[j]) {
		d2 = dist2(xgrid,ygrid,x[j],y[j],period);
		if(d2 < r2) {
		  /* point j covers grid point */
		  covered = TRUE;
		}
	      }
	    }
	  }
	  /* finished scanning all data points j */
	  if(!covered) kount = kount + 1;
	  /* finished considering grid point (xgrid,ygrid) */
	}
      }
    }
  }

  /*
    `kdisc' is the number of           grid points in the disc
    `kount' is the number of UNCOVERED grid points in the disc
  */

  if(areaint->hard) {
    if(kount == kdisc) cifval = 1.0;
    else cifval = 0.0;
  } else {
    /* usual calculation
       COVERED area fraction */
    covfrac = ((double) kdisc - (double) kount)/((double) kdisc);
    cifval = exp(areaint->logeta * covfrac);
  }

  return cifval;
}


Cifns AreaIntCifns = { &areaintInit, &areaintCif, (updafunptr) NULL, FALSE};
