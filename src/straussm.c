#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* for debugging code, include   #define DEBUG 1   */

/* Conditional intensity computation for Multitype Strauss process */

/* NOTE: types (marks) are numbered from 0 to ntypes-1 */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct MultiStrauss {
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

Cdata *straussminit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, ntypes, n2, m, hard;
  double g, r, r2, logg;
  MultiStrauss *multistrauss;

  multistrauss = (MultiStrauss *) R_alloc(1, sizeof(MultiStrauss));

  multistrauss->ntypes = ntypes = model.ntypes;
  n2 = ntypes * ntypes;

#ifdef DEBUG
  Rprintf("initialising space for %d types\n", ntypes);
#endif

  /* Allocate space for parameters */
  multistrauss->beta     = (double *) R_alloc((size_t) ntypes, sizeof(double));
  multistrauss->gamma    = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrauss->rad      = (double *) R_alloc((size_t) n2, sizeof(double));

  /* Allocate space for transformed parameters */
  multistrauss->rad2     = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrauss->loggamma = (double *) R_alloc((size_t) n2, sizeof(double));
  multistrauss->hard     = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Allocate scratch space for counts of each pair of types */
  multistrauss->kount     = (int *) R_alloc((size_t) n2, sizeof(int));

  /* Copy and process model parameters*/
  for(i = 0; i < ntypes; i++)
    multistrauss->beta[i]   = model.par[i];

  m = ntypes * (ntypes + 1);

  for(i = 0; i < ntypes; i++) {
    for(j = 0; j < ntypes; j++) {
      g = model.par[ntypes + i + j*ntypes];
      r = model.par[m + i + j*ntypes];
      r2 = r * r;
      hard = (g < DOUBLE_EPS);
      logg = (hard) ? 0 : log(g);
      MAT(multistrauss->gamma, i, j, ntypes) = g;
      MAT(multistrauss->rad, i, j, ntypes) = r;
      MAT(multistrauss->hard, i, j, ntypes) = hard; 
      MAT(multistrauss->loggamma, i, j, ntypes) = logg;
      MAT(multistrauss->rad2, i, j, ntypes) = r2;
    }
  }
  /* periodic boundary conditions? */
  multistrauss->period = model.period;
  multistrauss->per    = (model.period[0] > 0.0);

#ifdef DEBUG
  Rprintf("end initialiser\n");
#endif
  return((Cdata *) multistrauss);
}

/* conditional intensity evaluator */

double straussmcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ntypes, kount, ix, ixp1, j, mrk, mrkj, m1, m2;
  int *marks;
  double *x, *y;
  double u, v, lg;
  double d2, a, cifval;
  MultiStrauss *multistrauss;

  multistrauss = (MultiStrauss *) cdata;

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

  cifval = multistrauss->beta[mrk];

  if(npts == 0) 
    return(cifval);

  ntypes = multistrauss->ntypes;

#ifdef DEBUG
  Rprintf("initialising pair counts\n");
#endif

  /* initialise pair counts */
  for(m1 = 0; m1 < ntypes; m1++)
    for(m2 = 0; m2 < ntypes; m2++)
      MAT(multistrauss->kount, m1, m2, ntypes) = 0;

  /* compile pair counts */

#ifdef DEBUG
  Rprintf("compiling pair counts\n");
#endif

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(multistrauss->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],multistrauss->period);
	if(d2 < MAT(multistrauss->rad2, mrk, mrkj, ntypes)) 
	  MAT(multistrauss->kount, mrk, mrkj, ntypes)++;
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],multistrauss->period);
	if(d2 < MAT(multistrauss->rad2, mrk, mrkj, ntypes)) 
	  MAT(multistrauss->kount, mrk, mrkj, ntypes)++;
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	a = MAT(multistrauss->rad2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) MAT(multistrauss->kount, mrk, mrkj, ntypes)++;
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	a = MAT(multistrauss->rad2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j], 2);
	  if(a > 0) MAT(multistrauss->kount, mrk, mrkj, ntypes)++;
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
      kount = MAT(multistrauss->kount, m1, m2, ntypes);
      if(MAT(multistrauss->hard, m1, m2, ntypes)) {
	if(kount > 0) {
	  cifval = 0.0;
	  return(cifval);
	}
      } else {
	lg = MAT(multistrauss->loggamma, m1, m2, ntypes);
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
