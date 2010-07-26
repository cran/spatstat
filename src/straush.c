#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Hard core Strauss process */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct StraussHard {
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

Cdata *straushinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  StraussHard *strausshard;
  strausshard = (StraussHard *) R_alloc(1, sizeof(StraussHard));

  /* Interpret model parameters*/
  strausshard->beta   = model.par[0];
  strausshard->gamma  = model.par[1];
  strausshard->r      = model.par[2]; /* No longer passed as r^2 */
  strausshard->h      = model.par[3]; /* No longer passed as h^2 */
  strausshard->r2     = pow(strausshard->r, 2);
  strausshard->h2     = pow(strausshard->h, 2); 
  strausshard->r2h2   = strausshard->r2 - strausshard->h2;
  strausshard->period = model.period;
  /* is the interaction numerically equivalent to hard core ? */
  strausshard->hard   = (strausshard->gamma < DOUBLE_EPS);
  strausshard->loggamma = (strausshard->hard) ? 0 : log(strausshard->gamma);
  /* periodic boundary conditions? */
  strausshard->per    = (model.period[0] > 0.0);
  return((Cdata *) strausshard);
}

/* conditional intensity evaluator */

double straushcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, kount, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double r2, d2, h2, r2h2, a, cifval;
  StraussHard *strausshard;

  strausshard = (StraussHard *) cdata;

  period = strausshard->period;
  r2     = strausshard->r2;
  h2     = strausshard->h2;
  r2h2   = strausshard->r2h2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;

  npts = state.npts;

  if(npts == 0) 
    return(strausshard->beta);

  kount = 0;
  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(strausshard->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],strausshard->period);
	if(d2 < r2) {
	  if(d2 < h2) return(0.0);
	  kount = kount+1;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],strausshard->period);
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

  if(strausshard->hard) {
    if(kount > 0) cifval = 0.0;
    else cifval = strausshard->beta;
  }
  else cifval = strausshard->beta*exp(strausshard->loggamma*kount);
  
  return cifval;
}

Cifns StraussHardCifns = { &straushinit, &straushcif, (updafunptr) NULL, FALSE};
