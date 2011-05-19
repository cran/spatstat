#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Hard core process */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Hardcore {
  double beta;
  double h;   /* hard core distance */
  double h2;
  double *period;
  int per;
} Hardcore;


/* initialiser function */

Cdata *hardcoreinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Hardcore *hardcore;
  hardcore = (Hardcore *) R_alloc(1, sizeof(Hardcore));

  /* Interpret model parameters*/
  hardcore->beta   = model.par[0];
  hardcore->h      = model.par[1];
  hardcore->h2     = pow(hardcore->h, 2); 
  hardcore->period = model.period;
  /* periodic boundary conditions? */
  hardcore->per    = (model.period[0] > 0.0);

  return((Cdata *) hardcore);
}

/* conditional intensity evaluator */

double hardcorecif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, h2, a, cifval;
  Hardcore *hardcore;

  hardcore = (Hardcore *) cdata;

  h2     = hardcore->h2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return(hardcore->beta);

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(hardcore->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],hardcore->period);
	if(d2 < h2) return(0.0);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],hardcore->period);
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

  cifval = hardcore->beta;
  
  return cifval;
}

Cifns HardcoreCifns = { &hardcoreinit, &hardcorecif, (updafunptr) NULL, FALSE};
