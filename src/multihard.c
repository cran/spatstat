#include <R.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* for debugging code, include   #define DEBUG 1   */

/* Conditional intensity computation for Multitype Hardcore process */

/* NOTE: types (marks) are numbered from 0 to ntypes-1 */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct MultiHard {
  int ntypes;
  double *beta;    /* beta[i]  for i = 0 ... ntypes-1 */
  double *hc;      /* hc[i,j] = hc[j+ntypes*i] for i,j = 0... ntypes-1 */
  double *hc2;    /* squared radii */
  double *period;
  int per;
} MultiHard;


/* initialiser function */

Cdata *multihardinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j, ntypes, n2;
  double h;
  MultiHard *multihard;

  multihard = (MultiHard *) R_alloc(1, sizeof(MultiHard));

  multihard->ntypes = ntypes = model.ntypes;
  n2 = ntypes * ntypes;

#ifdef DEBUG
  Rprintf("initialising space for %d types\n", ntypes);
#endif

  /* Allocate space for parameters */
  multihard->beta     = (double *) R_alloc((size_t) ntypes, sizeof(double));
  multihard->hc       = (double *) R_alloc((size_t) n2, sizeof(double));

  /* Allocate space for transformed parameters */
  multihard->hc2      = (double *) R_alloc((size_t) n2, sizeof(double));

  /* Copy and process model parameters*/
  for(i = 0; i < ntypes; i++)
    multihard->beta[i]   = model.par[i];

  for(i = 0; i < ntypes; i++) {
    for(j = 0; j < ntypes; j++) {
      h = model.par[ntypes + i + j*ntypes];
      MAT(multihard->hc, i, j, ntypes) = h; 
      MAT(multihard->hc2, i, j, ntypes) = h * h;
    }
  }
  /* periodic boundary conditions? */
  multihard->period = model.period;
  multihard->per    = (model.period[0] > 0.0);

#ifdef DEBUG
  Rprintf("end initialiser\n");
#endif
  return((Cdata *) multihard);
}

/* conditional intensity evaluator */

double multihardcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ntypes, kount, ix, ixp1, j, mrk, mrkj, m1, m2;
  int *marks;
  double *x, *y;
  double u, v;
  double d2, a, cifval;
  MultiHard *multihard;

  multihard = (MultiHard *) cdata;

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

  cifval = multihard->beta[mrk];

  if(npts == 0) 
    return(cifval);

  ntypes = multihard->ntypes;

#ifdef DEBUG
  Rprintf("scanning data\n");
#endif

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(multihard->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],multihard->period);
	if(d2 < MAT(multihard->hc2, mrk, mrkj, ntypes)) {
	  cifval = 0;
	  return(cifval);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	d2 = dist2(u,v,x[j],y[j],multihard->period);
	if(d2 < MAT(multihard->hc2, mrk, mrkj, ntypes)) {
	  cifval = 0;
	  return(cifval);
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	mrkj = marks[j];
	a = MAT(multihard->hc2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j],2);
	  if(a > 0) {
	    cifval = 0;
	    return(cifval);
	  }
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	mrkj = marks[j];
	a = MAT(multihard->hc2, mrk, mrkj, ntypes); 
	a -= pow(u - x[j], 2);
	if(a > 0) {
	  a -= pow(v - y[j],2);
	  if(a > 0) {
	    cifval = 0;
	    return(cifval);
	  }
	}
      }
    }
  }

#ifdef DEBUG
  Rprintf("returning positive cif\n");
#endif
  return cifval;
}

Cifns MultiHardCifns = { &multihardinit, &multihardcif, (updafunptr) NULL, TRUE};
