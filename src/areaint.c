#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

#define NGRID 16

/*
  Conditional intensity function for an area-interaction process:

  cif = beta * eta^(1-B) where B = (uncovered area)/(pi r^2)

*/

/* Storage of parameters and precomputed/auxiliary data */

struct {
  /* model parameters */
  double beta;
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

void areaintInit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  double r, dx, dy, x0;
  int i, my, kdisc;
  /* Interpret model parameters*/
  AreaInt.beta   = model.par[0];
  AreaInt.eta    = model.par[1];
  AreaInt.r      = r = model.par[2]; 
  AreaInt.r2     = r * r;
  AreaInt.range2 = 4 * r * r;    /* square of interaction distance */
  /* is the model numerically equivalent to hard core ? */
  AreaInt.hard   = (AreaInt.eta < DOUBLE_EPS);
  AreaInt.logeta = (AreaInt.hard) ? 0 : log(AreaInt.eta);
  /* periodic boundary conditions? */
  AreaInt.period = model.period;
  AreaInt.per    = (model.period[0] > 0.0);
  /* grid counting */
  dx = dy = AreaInt.dx = (2 * r)/NGRID;
  AreaInt.xgrid0 = -r + dx/2;
  AreaInt.my = (int *) R_alloc((long) NGRID, sizeof(int));
  kdisc = 0;
  for(i = 0; i < NGRID; i++) {
    x0 = AreaInt.xgrid0 + i * dx;
    my = floor(sqrt(r * r - x0 * x0)/dy);
    my = (my < 0) ? 0 : my;
    AreaInt.my[i] = my;
    kdisc += 2 * my + 1;
  }
  AreaInt.kdisc = kdisc;
  /* allocate space for neighbour indicators */
  AreaInt.neigh = (char*) R_alloc((long) state.npmax, sizeof(char));
}

/* conditional intensity evaluator */

double areaintCif(prop, state)
     Propo prop;
     State state;
{
  int npts, kount, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double r2, d2, dx, dy, range2;
  double a, xgrid0, ygrid, ygrid0, xgrid, covfrac, cifval;
  int kdisc, kx, my, ky, covered;

  period = AreaInt.period;
  r2     = AreaInt.r2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;
  if(npts == 0) return AreaInt.beta;

  r2 = AreaInt.r2;
  dy = dx = AreaInt.dx;
  range2 = AreaInt.r2;    /* square of interaction distance */

  if(!AreaInt.per) {
    /*
      Euclidean distance
      First identify which data points are neighbours of (u,v)
    */
    ixp1 = ix + 1;
    /* If ix = NONE = -1, then ixp1 = 0 is correct */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	a = range2 - pow(u - x[j], 2);
	AreaInt.neigh[j] = FALSE;
	if(a > 0.) {
	  a -= pow(v - y[j], 2);
	  if(a > 0.)
	    AreaInt.neigh[j] = TRUE;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j < npts; j++) {
	a = range2 - pow(u - x[j], 2);
	AreaInt.neigh[j] = FALSE;
	if(a > 0.) {
	  a -= pow(v - y[j], 2);
	  if(a > 0.)
	    AreaInt.neigh[j] = TRUE;
	}
      }
    }
    /* scan a grid of points centred at (u,v) */
    kount = 0;
    xgrid0 = u + AreaInt.xgrid0;
    for(kx=0; kx<NGRID; kx++) {
      xgrid = xgrid0 + kx * dx;
      my = AreaInt.my[kx];
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
	    if(AreaInt.neigh[j]) {
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
	    if(AreaInt.neigh[j]) {
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
    kdisc = AreaInt.kdisc;
  } else {
    /*
      periodic distance
      First identify which data points are neighbours of (u,v)
    */
    ixp1 = ix + 1;
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],period);
	AreaInt.neigh[j] = (d2 < range2);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],period);
	AreaInt.neigh[j] = (d2 < range2);
      }
    }
    /* scan a grid of ngrid * ngrid points centred at (u,v) */
    kount = 0;
    kdisc = 0;
    xgrid0 = u + AreaInt.xgrid0;
    ygrid0 = u + AreaInt.xgrid0;
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
	      if(AreaInt.neigh[j]) {
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
	      if(AreaInt.neigh[j]) {
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

  if(AreaInt.hard) {
    if(kount == kdisc) cifval = AreaInt.beta;
    else cifval = 0.0;
  } else {
    /* usual calculation
       COVERED area fraction */
    covfrac = ((double) kdisc - (double) kount)/((double) kdisc);
    cifval = AreaInt.beta * exp(AreaInt.logeta * covfrac);
  }

  return cifval;
}


Cifns AreaIntCifns = { &areaintInit, &areaintCif, (updafunptr) NULL, FALSE};
