#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Hard core process */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  double beta;
  double h;   /* hard core distance */
  double h2;
  double *period;
  int per;
} Hardcore;


/* initialiser function */

void hardcoreinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  Hardcore.beta   = model.par[0];
  Hardcore.h      = model.par[1];
  Hardcore.h2     = pow(Hardcore.h, 2); 
  Hardcore.period = model.period;
  /* periodic boundary conditions? */
  Hardcore.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double hardcorecif(prop, state)
     Propo prop;
     State state;
{
  int npts, kount, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double d2, h2, a, cifval;

  period = Hardcore.period;
  h2     = Hardcore.h2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return(Hardcore.beta);

  kount = 0;
  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(Hardcore.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],Hardcore.period);
	if(d2 < h2) return(0.0);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Hardcore.period);
	if(d2 < h2) return(0.0);
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	a = h2 - pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) 
	    return(0.0);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	a = h2 - pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0)
	    return(0.0);
	}
      }
    }
  }

  cifval = Hardcore.beta;
  
  return cifval;
}

Cifns HardcoreCifns = { &hardcoreinit, &hardcorecif, (updafunptr) NULL, FALSE};
