#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* for debugging code, include   #define DEBUG 1   */

/* Conditional intensity computation for Multitype Strauss process */

/* NOTE: types (marks) are numbered from 0 to ntypes-1 */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  int ntypes;
  double *beta;    /* beta[i]  for i = 0 ... ntypes-1 */
  double *gamma;   /* gamma[i,j] = gamma[i+ntypes*j] for i,j = 0... ntypes-1 */
  double *rad;     /* rad[i,j] = rad[j+ntypes*i] for i,j = 0... ntypes-1 */
  double *rad2;    /* squared radii */
  double *loggamma; /* logs of gamma[i,j] */
  double *period;
  int    *hard;     /* hard[i,j] = 1 if gamma[i,j] ~~ 0 */
  int    *kount;    /* space for kounting pairs of each type */
  int per;
} MultiStrauss;


/* initialiser function */

void straussminit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, ntypes, n2, m, hard;
  double g, r, r2, logg;

  MultiStrauss.ntypes = ntypes = model.ntypes;
  n2 = ntypes * ntypes;

#ifdef DEBUG
  Rprintf("initialising space for %d types\n", ntypes);
#endif

  /* Allocate space for parameters */
  MultiStrauss.beta     = (double *) R_alloc((size_t) ntypes, sizeof(double));
  MultiStrauss.gamma    = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStrauss.rad      = (double *) R_alloc((size_t) n2, sizeof(double));

  /* Allocate space for transformed parameters */
  MultiStrauss.rad2     = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStrauss.loggamma = (double *) R_alloc((size_t) n2, sizeof(double));
  MultiStrauss.hard     = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Allocate scratch space for counts of each pair of types */
  MultiStrauss.kount     = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Copy and process model parameters*/
  for(i = 0; i < ntypes; i++)
    MultiStrauss.beta[i]   = model.par[i];

  m = ntypes * (ntypes + 1);

  for(i = 0; i < ntypes; i++) {
    for(j = 0; j < ntypes; j++) {
      g = model.par[ntypes + i + j*ntypes];
      r = model.par[m + i + j*ntypes];
      r2 = r * r;
      hard = (g < DOUBLE_EPS);
      logg = (hard) ? 0 : log(g);
      MAT(MultiStrauss.gamma, i, j, ntypes) = g;
      MAT(MultiStrauss.rad, i, j, ntypes) = r;
      MAT(MultiStrauss.hard, i, j, ntypes) = hard; 
      MAT(MultiStrauss.loggamma, i, j, ntypes) = logg;
      MAT(MultiStrauss.rad2, i, j, ntypes) = r2;
    }
  }
  /* periodic boundary conditions? */
  MultiStrauss.period = model.period;
  MultiStrauss.per    = (model.period[0] > 0.0);

#ifdef DEBUG
  Rprintf("end initialiser\n");
#endif
}

/* conditional intensity evaluator */

double straussmcif(prop, state)
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

  cifval = MultiStrauss.beta[mrk];

  if(npts == 0) 
    return(cifval);

  ntypes = MultiStrauss.ntypes;

#ifdef DEBUG
  Rprintf("initialising pair counts\n");
#endif

  /* initialise pair counts */
  for(m1 = 0; m1 < ntypes; m1++)
    for(m2 = 0; m2 < ntypes; m2++)
      MAT(MultiStrauss.kount, m1, m2, ntypes) = 0;

  /* compile pair counts */

#ifdef DEBUG
  Rprintf("compiling pair counts\n");
#endif

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(MultiStrauss.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],MultiStrauss.period);
	if(d2 < MAT(MultiStrauss.rad2, mrk, mrkj, ntypes)) 
	  MAT(MultiStrauss.kount, mrk, mrkj, ntypes)++;
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],MultiStrauss.period);
	if(d2 < MAT(MultiStrauss.rad2, mrk, mrkj, ntypes)) 
	  MAT(MultiStrauss.kount, mrk, mrkj, ntypes)++;
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	a = MAT(MultiStrauss.rad2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) MAT(MultiStrauss.kount, mrk, mrkj, ntypes)++;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	a = MAT(MultiStrauss.rad2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) MAT(MultiStrauss.kount, mrk, mrkj, ntypes)++;
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
      kount = MAT(MultiStrauss.kount, m1, m2, ntypes);
      if(MAT(MultiStrauss.hard, m1, m2, ntypes)) {
	if(kount > 0) {
	  cifval = 0.0;
	  return(cifval);
	}
      } else {
	lg = MAT(MultiStrauss.loggamma, m1, m2, ntypes);
	cifval *= exp(lg * kount);
      }
    }
  }
  
#ifdef DEBUG
  Rprintf("returning positive cif\n");
#endif
  return cifval;
}

Cifns MultiStraussCifns = { &straussminit, &straussmcif, (updafunptr) NULL, TRUE};
