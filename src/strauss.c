#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Strauss process */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  double beta;
  double gamma;
  double r;
  double loggamma;
  double r2;
  double *period;
  int hard;
  int per;
} Strauss;


/* initialiser function */

void straussinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  Strauss.beta   = model.par[0];
  Strauss.gamma  = model.par[1];
  Strauss.r      = model.par[2]; /* No longer passed as r^2 */
  Strauss.r2     = Strauss.r * Strauss.r; 
  Strauss.period = model.period;
  /* is the model numerically equivalent to hard core ? */
  Strauss.hard   = (Strauss.gamma < DOUBLE_EPS);
  Strauss.loggamma = (Strauss.hard) ? 0 : log(Strauss.gamma);
  /* periodic boundary conditions? */
  Strauss.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double strausscif(prop, state)
     Propo prop;
     State state;
{
  int npts, kount, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double r2, d2, a, cifval;

  period = Strauss.period;
  r2     = Strauss.r2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return(Strauss.beta);

  kount = 0;
  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(Strauss.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],Strauss.period);
	if(d2 < r2) kount = kount+1;
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Strauss.period);
	if(d2 < r2) kount = kount+1;
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	a = r2 - pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) kount = kount+1;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	a = r2 - pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) kount = kount+1;
	}
      }
    }
  }

  if(Strauss.hard) {
    if(kount > 0) cifval = 0.0;
    else cifval = Strauss.beta;
  }
  else cifval = Strauss.beta*exp(Strauss.loggamma*kount);
  
  return cifval;
}

Cifns StraussCifns = { &straussinit, &strausscif, (updafunptr) NULL, FALSE};
