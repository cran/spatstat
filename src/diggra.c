#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"

/* Conditional intensity computation for Diggle-Gratton process */

/*
 Conditional intensity function for a pairwise interaction point
 process with interaction function as given by 

                  e(t) = 0 for t < delta
                       = (t-delta)/(rho-delta)^kappa for delta <= t < rho
                       = 1 for t >= rho

 (See page 767 of Diggle, Gates, and Stibbard, Biometrika vol. 74,
  1987, pages 763 -- 770.)
*/

/* Storage of parameters and precomputed/auxiliary data */

typedef struct Diggra {
  double beta;
  double kappa;
  double delta;
  double rho;
  double delta2;  /*  delta^2   */
  double rho2;    /*  rho^2 */
  double fac;   /*   1/(rho-delta)  */
  double *period;
  int per;
} Diggra;


/* initialiser function */

Cdata *diggrainit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Diggra *diggra;
  diggra = (Diggra *) R_alloc(1, sizeof(Diggra));

  /* Interpret model parameters*/
  diggra->beta   = model.par[0];
  diggra->kappa  = model.par[1];
  diggra->delta  = model.par[2];
  diggra->rho    = model.par[3];
  diggra->period = model.period;
  /* constants */
  diggra->delta2 = pow(diggra->delta, 2);
  diggra->rho2 = pow(diggra->rho, 2);
  diggra->fac = 1/(diggra->rho - diggra->delta);
  /* periodic boundary conditions? */
  diggra->per    = (model.period[0] > 0.0);
  return((Cdata *) diggra);
}

/* conditional intensity evaluator */

double diggracif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *x, *y;
  double u, v;
  double d2, pairprod, cifval;
  Diggra *diggra;

  diggra = (Diggra *) cdata;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = diggra->beta;

  if(npts == 0) 
    return(cifval);

  pairprod = 1;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(diggra->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],diggra->period);
	if(d2 < diggra->delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < diggra->rho2) {
	  pairprod *= diggra->fac * (sqrt(d2)-diggra->delta);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],diggra->period);
	if(d2 < diggra->delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < diggra->rho2) {
	  pairprod *= diggra->fac * (sqrt(d2)-diggra->delta);
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < diggra->delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < diggra->rho2) {
	  pairprod *= diggra->fac * (sqrt(d2)-diggra->delta);
	}
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < diggra->delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < diggra->rho2) {
	  pairprod *= diggra->fac * (sqrt(d2)-diggra->delta);
	}
      }
    }
  }

  cifval *= pow(pairprod, diggra->kappa);
  return cifval;
}

Cifns DiggraCifns = { &diggrainit, &diggracif, (updafunptr) NULL, FALSE};

