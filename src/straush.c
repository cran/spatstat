#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Hard core Strauss process */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  double beta;
  double gamma;
  double r;   /* interaction distance */
  double h;   /* hard core distance */
  double loggamma;
  double r2;
  double h2;
  double r2h2;  /* r^2 - h^2 */
  double *period;
  int hard;
  int per;
} StraussHard;


/* initialiser function */

void straushinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  StraussHard.beta   = model.par[0];
  StraussHard.gamma  = model.par[1];
  StraussHard.r      = model.par[2]; /* No longer passed as r^2 */
  StraussHard.h      = model.par[3]; /* No longer passed as h^2 */
  StraussHard.r2     = pow(StraussHard.r, 2);
  StraussHard.h2     = pow(StraussHard.h, 2); 
  StraussHard.r2h2   = StraussHard.r2 - StraussHard.h2;
  StraussHard.period = model.period;
  /* is the interaction numerically equivalent to hard core ? */
  StraussHard.hard   = (StraussHard.gamma < DOUBLE_EPS);
  StraussHard.loggamma = (StraussHard.hard) ? 0 : log(StraussHard.gamma);
  /* periodic boundary conditions? */
  StraussHard.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double straushcif(prop, state)
     Propo prop;
     State state;
{
  int npts, kount, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double r2, d2, h2, r2h2, a, cifval;

  period = StraussHard.period;
  r2     = StraussHard.r2;
  h2     = StraussHard.h2;
  r2h2   = StraussHard.r2h2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return(StraussHard.beta);

  kount = 0;
  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(StraussHard.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],StraussHard.period);
	if(d2 < r2) {
	  if(d2 < h2) return(0.0);
	  kount = kount+1;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],StraussHard.period);
	if(d2 < r2) {
	  if(d2 < h2) return(0.0);
	  kount = kount+1;
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	a = r2 - pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) {
	    if(a > r2h2)
	      return(0.0);
	    kount = kount+1;
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	a = r2 - pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) {
	    if(a > r2h2)
	      return(0.0);
	    kount = kount+1;
	  }
	}
      }
    }
  }

  if(StraussHard.hard) {
    if(kount > 0) cifval = 0.0;
    else cifval = StraussHard.beta;
  }
  else cifval = StraussHard.beta*exp(StraussHard.loggamma*kount);
  
  return cifval;
}

Cifns StraussHardCifns = { &straushinit, &straushcif, (updafunptr) NULL, FALSE};
