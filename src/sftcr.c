#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"


/* Conditional intensity computation for Soft Core process */

/* Storage of parameters and precomputed/auxiliary data */

struct {
  double beta;
  double sigma;
  double kappa;
  double ook;  /*   1/kappa     */
  double stuk; /* sigma^(2/kappa) */
  double *period;
  int per;
} Softcore;


/* initialiser function */

void sftcrinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  Softcore.beta   = model.par[0];
  Softcore.sigma  = model.par[1];
  Softcore.kappa  = model.par[2];
  Softcore.period = model.period;
  /* constants */
  Softcore.ook = 1/Softcore.kappa;
  Softcore.stuk = pow(Softcore.sigma, 2/Softcore.kappa);
  /* periodic boundary conditions? */
  Softcore.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double sftcrcif(prop, state)
     Propo prop;
     State state;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, pairsum, cifval, ook, stuk;

  ook = Softcore.ook;
  stuk = Softcore.stuk;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = Softcore.beta;

  if(npts == 0) 
    return(cifval);

  pairsum = 0;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(Softcore.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],Softcore.period);
	pairsum += pow(d2, ook);
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Softcore.period);
	pairsum += pow(d2, ook);
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	pairsum += pow(d2, ook);
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j],2) + pow(v-y[j],2);
	pairsum += pow(d2, ook);
      }
    }
  }

  cifval *= exp(-stuk * pairsum);
  return cifval;
}

Cifns SoftcoreCifns = { &sftcrinit, &sftcrcif, (updafunptr) NULL, FALSE};

