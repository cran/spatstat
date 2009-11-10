#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* for debugging code, include   #define DEBUG 1   */

/* Conditional intensity computation for Multitype Strauss hardcore process */

/* NOTE: types (marks) are numbered from 0 to ntypes-1 */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  int ntypes;
  double *beta;    /* beta[i]  for i = 0 ... ntypes-1 */
  double *gamma;   /* gamma[i,j] = gamma[i+ntypes*j] for i,j = 0... ntypes-1 */
  double *rad;     /* rad[i,j] = rad[j+ntypes*i] for i,j = 0... ntypes-1 */
  double *hc;      /* hc[i,j] = hc[j+ntypes*i] for i,j = 0... ntypes-1 */
  double *rad2;    /* squared radii */
  double *hc2;    /* squared radii */
  double *rad2hc2;    /* r^2 - h^2 */
  double *loggamma; /* logs of gamma[i,j] */
  double *period;
  int    *hard;     /* hard[i,j] = 1 if gamma[i,j] ~~ 0 */
  int    *kount;    /* space for kounting pairs of each type */
  int per;
} MultiStraussHard;


/* initialiser function */

void straushminit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, ntypes, n2, m, mm, hard;
  double g, r, h, r2, h2, logg;

  MultiStraussHard.ntypes = ntypes = model.ntypes;
  n2 = ntypes * ntypes;

#ifdef DEBUG
  Rprintf("initialising space for %d types\n", ntypes);
#endif

  /* Allocate space for parameters */
  MultiStraussHard.beta     = (double *) R_alloc((size_t) ntypes, sizeof(double));
  MultiStraussHard.gamma    = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStraussHard.rad      = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStraussHard.hc       = (double *) R_alloc((size_t) n2, sizeof(double));

  /* Allocate space for transformed parameters */
  MultiStraussHard.rad2     = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStraussHard.hc2      = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStraussHard.rad2hc2  = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStraussHard.loggamma = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStraussHard.hard     = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Allocate scratch space for counts of each pair of types */
  MultiStraussHard.kount    = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Copy and process model parameters*/
  for(i = 0; i < ntypes; i++)
    MultiStraussHard.beta[i]   = model.par[i];

  m = ntypes * (ntypes + 1);
  mm = m + ntypes * ntypes;

  for(i = 0; i < ntypes; i++) {
    for(j = 0; j < ntypes; j++) {
      g = model.par[ntypes + i + j*ntypes];
      r = model.par[m + i + j*ntypes];
      h = model.par[mm + i + j*ntypes];
      r2 = r * r;
      h2 = h * h;
      hard = (g < DOUBLE_EPS);
      logg = (hard) ? 0 : log(g);
      MAT(MultiStraussHard.gamma, i, j, ntypes) = g;
      MAT(MultiStraussHard.rad, i, j, ntypes) = r;
      MAT(MultiStraussHard.hc, i, j, ntypes) = h; 
      MAT(MultiStraussHard.rad2, i, j, ntypes) = r2;
      MAT(MultiStraussHard.hc2, i, j, ntypes) = h2;
      MAT(MultiStraussHard.rad2hc2, i, j, ntypes) = r2-h2;
      MAT(MultiStraussHard.hard, i, j, ntypes) = hard; 
      MAT(MultiStraussHard.loggamma, i, j, ntypes) = logg;
    }
  }
  /* periodic boundary conditions? */
  MultiStraussHard.period = model.period;
  MultiStraussHard.per    = (model.period[0] > 0.0);

#ifdef DEBUG
  Rprintf("end initialiser\n");
#endif
}

/* conditional intensity evaluator */

double straushmcif(prop, state)
     Propo prop;
     State state;
{
  int npts, ntypes, kount, ix, ixp1, j, mrk, mrkj, m1, m2;
  int *marks;
  double *x, *y;
  double u, v, lg;
  double d2, a, cifval;

  u  = prop.u;
  v  = prop.v;
  mrk = prop.mrk;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  marks = state.marks;

  npts = state.npts;

#ifdef DEBUG
  Rprintf("computing cif: u=%lf, v=%lf, mrk=%d\n", u, v, mrk);
#endif

  cifval = MultiStraussHard.beta[mrk];

  if(npts == 0) 
    return(cifval);

  ntypes = MultiStraussHard.ntypes;

#ifdef DEBUG
  Rprintf("initialising pair counts\n");
#endif

  /* initialise pair counts */
  for(m1 = 0; m1 < ntypes; m1++)
    for(m2 = 0; m2 < ntypes; m2++)
      MAT(MultiStraussHard.kount, m1, m2, ntypes) = 0;

  /* compile pair counts */

#ifdef DEBUG
  Rprintf("compiling pair counts\n");
#endif

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(MultiStraussHard.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],MultiStraussHard.period);
	if(d2 < MAT(MultiStraussHard.rad2, mrk, mrkj, ntypes)) {
	  if(d2 < MAT(MultiStraussHard.hc2, mrk, mrkj, ntypes)) {
	    cifval = 0;
	    return(cifval);
	  }
	  MAT(MultiStraussHard.kount, mrk, mrkj, ntypes)++;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],MultiStraussHard.period);
	if(d2 < MAT(MultiStraussHard.rad2, mrk, mrkj, ntypes)) {
	  if(d2 < MAT(MultiStraussHard.hc2, mrk, mrkj, ntypes)) {
	    cifval = 0;
	    return(cifval);
	  }
	  MAT(MultiStraussHard.kount, mrk, mrkj, ntypes)++;
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	a = MAT(MultiStraussHard.rad2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j],2);
	  if(a > 0) {
	    if(a > MAT(MultiStraussHard.rad2hc2, mrk, mrkj, ntypes)) {
	      cifval = 0;
	      return(cifval);
	    }
	    MAT(MultiStraussHard.kount, mrk, mrkj, ntypes)++;
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	a = MAT(MultiStraussHard.rad2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j],2);
	  if(a > 0) {
	    if(a > MAT(MultiStraussHard.rad2hc2, mrk, mrkj, ntypes)) {
	      cifval = 0;
	      return(cifval);
	    }
	    MAT(MultiStraussHard.kount, mrk, mrkj, ntypes)++;
	  }
	}
      }
    }
  }

#ifdef DEBUG
  Rprintf("multiplying cif factors\n");
#endif
  /* multiply cif value by pair potential */
  for(m1 = 0; m1 < ntypes; m1++) {
    for(m2 = 0; m2 < ntypes; m2++) {
      kount = MAT(MultiStraussHard.kount, m1, m2, ntypes);
      if(MAT(MultiStraussHard.hard, m1, m2, ntypes)) {
	if(kount > 0) {
	  cifval = 0.0;
	  return(cifval);
	}
      } else {
	lg = MAT(MultiStraussHard.loggamma, m1, m2, ntypes);
	cifval *= exp(lg * kount);
      }
    }
  }
  
#ifdef DEBUG
  Rprintf("returning positive cif\n");
#endif
  return cifval;
}

Cifns MultiStraussHardCifns = { &straushminit, &straushmcif, (updafunptr) NULL, TRUE};
