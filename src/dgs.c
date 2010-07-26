#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "methas.h"
#include "dist2.h"


/* Conditional intensity computation for Diggle-Gates-Stibbard process */

/*
 Conditional intensity function for a pairwise interaction point
 process with interaction function as given by 

                  e(t) = sin^2(pi*t/2*rho) for t < rho
                       = 1 for t >= rho

 (See page 767 of Diggle, Gates, and Stibbard, Biometrika vol. 74,
  1987, pages 763 -- 770.)
*/

#define PION2 M_PI_2   /* pi/2 defined in Rmath.h */


/* Storage of parameters and precomputed/auxiliary data */

typedef struct Dgs {
  double beta;
  double rho;
  double rho2;
  double pion2rho;
  double *period;
  int per;
} Dgs;


/* initialiser function */

Cdata *dgsinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  Dgs *dgs;
  /* allocate storage */
  dgs = (Dgs *) R_alloc(1, sizeof(Dgs));

  /* Interpret model parameters*/
  dgs->beta   = model.par[0];
  dgs->rho    = model.par[1];
  dgs->period = model.period;
  /* constants */
  dgs->rho2       = pow(dgs->rho, 2);
  dgs->pion2rho   = PION2/dgs->rho;
  /* periodic boundary conditions? */
  dgs->per    = (model.period[0] > 0.0);
  return((Cdata *) dgs);
}

/* conditional intensity evaluator */

double dgscif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int npts, ix, ixp1, j;
  double *period, *x, *y;
  double u, v;
  double d2, r2, pairprod, cifval;
  Dgs *dgs;

  dgs = (Dgs *) cdata;

  period = dgs->period;
  r2 = dgs->rho2;

  u  = prop.u;
  v  = prop.v;
  ix = prop.ix;
  x  = state.x;
  y  = state.y;
  npts = state.npts;

  cifval = dgs->beta;

  if(npts == 0) 
    return(cifval);

  pairprod = 1;

  ixp1 = ix+1;
  /* If ix = NONE = -1, then ixp1 = 0 is correct */
  if(dgs->per) { /* periodic distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = dist2(u,v,x[j],y[j],dgs->period);
	if(d2 < r2) pairprod *= sin(dgs->pion2rho * sqrt(d2));
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],dgs->period);
	if(d2 < r2) pairprod *= sin(dgs->pion2rho * sqrt(d2));
      }
    }
  }
  else { /* Euclidean distance */
    if(ix > 0) {
      for(j=0; j < ix; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < r2) pairprod *= sin(dgs->pion2rho * sqrt(d2));
      }
    }
    if(ixp1 < npts) {
      for(j=ixp1; j<npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < r2) pairprod *= sin(dgs->pion2rho * sqrt(d2));
      }
    }
  }

  cifval *= pairprod * pairprod;
  return cifval;
}

Cifns DgsCifns = { &dgsinit, &dgscif, (updafunptr) NULL, FALSE};

