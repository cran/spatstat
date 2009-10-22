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

struct {
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

void diggrainit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  /* Interpret model parameters*/
  Diggra.beta   = model.par[0];
  Diggra.kappa  = model.par[1];
  Diggra.delta  = model.par[2];
  Diggra.rho    = model.par[3];
  Diggra.period = model.period;
  /* constants */
  Diggra.delta2 = pow(Diggra.delta, 2);
  Diggra.rho2 = pow(Diggra.rho, 2);
  Diggra.fac = 1/(Diggra.rho - Diggra.delta);
  /* periodic boundary conditions? */
  Diggra.per    = (model.period[0] > 0.0);
}

/* conditional intensity evaluator */

double diggracif(prop, state)
     Propo prop;
     State state;
{
  int npts, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double d2, r2, a, pairprod, cifval;

  period = Diggra.period;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = Diggra.beta;

  if(npts == 0) 
    return(cifval);

  pairprod = 1;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(Diggra.per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],Diggra.period);
	if(d2 < Diggra.delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Diggra.rho2) {
	  pairprod *= Diggra.fac * (sqrt(d2)-Diggra.delta);
	}
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Diggra.period);
	if(d2 < Diggra.delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Diggra.rho2) {
	  pairprod *= Diggra.fac * (sqrt(d2)-Diggra.delta);
	}
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < Diggra.delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Diggra.rho2) {
	  pairprod *= Diggra.fac * (sqrt(d2)-Diggra.delta);
	}
      }
    }  
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < Diggra.delta2) {
	  cifval = 0;
	  return(cifval);
	} else if(d2 < Diggra.rho2) {
	  pairprod *= Diggra.fac * (sqrt(d2)-Diggra.delta);
	}
      }
    }
  }

  cifval *= pow(pairprod, Diggra.kappa);
  return cifval;
}

Cifns DiggraCifns = { &diggrainit, &diggracif, (updafunptr) NULL, FALSE};

