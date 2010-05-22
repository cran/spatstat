#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Fiksel process */

/*
 Conditional intensity function for a pairwise interaction point
 process with interaction function 

                  e(t) = 0 for t < h
                       = exp(a * exp(- kappa * t)) for h <= t < r
                       = 1 for t >= r

*/

/* Storage of parameters and precomputed/auxiliary data */

struct {
  double beta;
  double r;
  double h;
  double kappa;
  double a;
  double h2;  /*  h^2   */
  double r2;  /*  r^2 */
  double *period;
  int per;
} Fiksel;


/* initialiser function */

void fikselinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  Fiksel.beta   = model.par[0];
  Fiksel.r      = model.par[1];
  Fiksel.h      = model.par[2];
  Fiksel.kappa  = model.par[3];
  Fiksel.a      = model.par[4];
  Fiksel.period = model.period;
  /* constants */
  Fiksel.h2 = pow(Fiksel.h, 2);
  Fiksel.r2 = pow(Fiksel.r, 2);
  /* periodic boundary conditions? */
  Fiksel.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double fikselcif(prop, state)
     Propo prop;
     State state;
{
  int npts, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double d2, pairpotsum, cifval;

  period = Fiksel.period;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = 0;

  if(npts == 0) 
    return(cifval);

  pairpotsum = 0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(Fiksel.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],Fiksel.period);
	if(d2 < Fiksel.h2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Fiksel.r2) {
	  pairpotsum += exp(-Fiksel.kappa * sqrt(d2));
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Fiksel.period);
	if(d2 < Fiksel.h2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Fiksel.r2) {
	  pairpotsum += exp(-Fiksel.kappa * sqrt(d2));
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < Fiksel.h2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Fiksel.r2) {
	  pairpotsum += exp(-Fiksel.kappa * sqrt(d2));
	}
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < Fiksel.h2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Fiksel.r2) {
	  pairpotsum += exp(-Fiksel.kappa * sqrt(d2));
	}
      }
    }
  }

  cifval = Fiksel.beta * exp(Fiksel.a * pairpotsum);
  return cifval;
}

Cifns FikselCifns = { &fikselinit, &fikselcif, (updafunptr) NULL, FALSE};

