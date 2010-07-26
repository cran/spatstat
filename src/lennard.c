#include <R.h>
#include <Rmath.h>
#include <R_ext/Constants.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Lennard-Jones process */

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Lennard {
  double beta;
  double sigma;
  double epsilon;
  double sigma2;  /*   sigma^2     */
  double foureps;    /*   4 * epsilon     */
  double d2min;  /* minimum value of d^2 which yields nonzero intensity */
  double d2max;  /* maximum value of d^2 which has nontrivial contribution */
  double *period;
  int per;
} Lennard;

/* 
   This is intended to be the largest x such that exp(-x) != 0 
   although the exact value is not needed
*/
#define MAXEXP (-log(DOUBLE_XMIN))
#define MINEXP (log(1.001))

/* initialiser function */

Cdata *lennardinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Lennard *lennard;

  lennard = (Lennard *) R_alloc(1, sizeof(Lennard));

  /* Interpret model parameters*/
  lennard->beta    = model.par[0];
  lennard->sigma   = model.par[1];
  lennard->epsilon = model.par[2];
  lennard->period  = model.period;
  /* constants */
  lennard->sigma2  = pow(lennard->sigma, 2);
  lennard->foureps = 4 * lennard->epsilon;
  lennard->d2min   = lennard->sigma2 * pow(lennard->foureps/MAXEXP, 1/6);
  lennard->d2max   = lennard->sigma2 * pow(lennard->foureps/MINEXP, 1/3);
  /* periodic boundary conditions? */
  lennard->per    = (model.period[0] > 0.0);

  return((Cdata *) lennard);
}

/* conditional intensity evaluator */

double lennardcif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, sigma2, ratio6, pairsum, cifval;
  Lennard *lennard;

  lennard = (Lennard *) cdata;

  sigma2 = lennard->sigma2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = lennard->beta;

  if(npts == 0) 
    return(cifval);

  pairsum = 0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(lennard->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],lennard->period);
	if(d2 < lennard->d2max) {
	  if(d2 < lennard->d2min) {
	    cifval = 0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 - pow(ratio6, 2);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],lennard->period);
	if(d2 < lennard->d2max) {
	  if(d2 < lennard->d2min) {
	    cifval = 0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 - pow(ratio6, 2);
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	if(d2 < lennard->d2max) {
	  if(d2 < lennard->d2min) {
	    cifval = 0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 - pow(ratio6, 2);
	}
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	if(d2 < lennard->d2max) {
	  if(d2 < lennard->d2min) {
	    cifval = 0;
	    return cifval;
	  }
	  ratio6 = pow(sigma2/d2, 3);
	  pairsum += ratio6 - pow(ratio6, 2);
	}
      }
    }
  }

  cifval *= exp(lennard->foureps * pairsum);
  return cifval;
}

Cifns LennardCifns = { &lennardinit, &lennardcif, (updafunptr) NULL, FALSE};

