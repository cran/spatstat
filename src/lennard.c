#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"


/* Conditional intensity computation for Lennard-Jones process */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  double beta;
  double sigma;
  double epsilon;
  double sigma2;  /*   sigma^2     */
  double foureps;    /*   4 * epsilon     */
  double *period;
  int per;
} Lennard;


/* initialiser function */

void lennardinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  Lennard.beta    = model.par[0];
  Lennard.sigma   = model.par[1];
  Lennard.epsilon = model.par[2];
  Lennard.period  = model.period;
  /* constants */
  Lennard.foureps = 4 * Lennard.epsilon;
  /* periodic boundary conditions? */
  Lennard.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double lennardcif(prop, state)
     Propo prop;
     State state;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, sigma2, ratio6, pairsum, cifval, foureps;

  foureps = Lennard.foureps;
  sigma2 = Lennard.sigma2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = Lennard.beta;

  if(npts == 0) 
    return(cifval);

  pairsum = 0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(Lennard.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],Lennard.period);
	ratio6 = pow(sigma2/d2, 3);
	pairsum += ratio6 - pow(ratio6, 2);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Lennard.period);
	ratio6 = pow(sigma2/d2, 3);
	pairsum += ratio6 - pow(ratio6, 2);
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	ratio6 = pow(sigma2/d2, 3);
	pairsum += ratio6 - pow(ratio6, 2);
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	ratio6 = pow(sigma2/d2, 3);
	pairsum += ratio6 - pow(ratio6, 2);
      }
    }
  }

  cifval *= exp(foureps * pairsum);
  return cifval;
}

Cifns LennardCifns = { &lennardinit, &lennardcif, (updafunptr) NULL, FALSE};

