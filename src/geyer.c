#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

void fexitc(const char *msg);

/*
  Conditional intensity function for a Geyer saturation process.  
*/

struct {
  /* model parameters */
  double beta;
  double gamma;
  double r;
  double s;
  /* transformations of the parameters */
  double r2;
  double loggamma;
  int hard;
  /* periodic distance */
  double *period;
  int per;
  /* auxiliary counts */
  int *aux;
} Geyer;

void geyerinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j;
  double r, d2;
  /* Interpret model parameters*/
  Geyer.beta   = model.par[0];
  Geyer.gamma  = model.par[1];
  Geyer.r      = model.par[2]; /* not squared any more */
  Geyer.s      = model.par[3]; 
  Geyer.r2     = Geyer.r * Geyer.r;
  /* is the model numerically equivalent to hard core ? */
  Geyer.hard   = (Geyer.gamma < DOUBLE_EPS);
  Geyer.loggamma = (Geyer.hard) ? 0 : log(Geyer.gamma);
  /* periodic boundary conditions? */
  Geyer.period = model.period;
  Geyer.per    = (model.period[0] > 0.0);
  /* allocate storage for auxiliary counts */
  Geyer.aux = (int *) R_alloc((size_t) state.npmax, sizeof(int));
  /* Initialise auxiliary counts */
  for(i = 0; i < state.npmax; i++) 
    Geyer.aux[i] = 0;

  for(i = 0; i < state.npts; i++) {
    for(j = 0; j < state.npts; j++) {
      d2 = dist2either(state.x[i], state.y[i], state.x[j], state.y[j], 
		       Geyer.period);
      if(d2 < Geyer.r2)
	Geyer.aux[i] += 1;
    }
  }
}

double geyercif(prop, state)
     Propo prop;
     State state;
{
  int ix, j, npts, tee;
  double u, v, d2, r2, s;
  double w, a, dd2, b, f, cifval;
  double *x, *y;
  int *aux;

  npts = state.npts;
  if(npts==0) return Geyer.beta;

  x = state.x;
  y = state.y;
  u = prop.u;
  v = prop.v;
  ix = prop.ix;

  r2 = Geyer.r2;
  s = Geyer.s;

  aux = Geyer.aux;

  /* 
     tee = neighbour count at the point in question;
     w   = sum of changes in (saturated) neighbour counts at other points 
  */
  tee = w = 0.0;

  if(prop.itype == BIRTH) {
    if(Geyer.per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Geyer.period);
	if(d2 < r2) {
	  tee++;
	  f = s - aux[j];
	  if(f > 1) /* j is not saturated after addition of (u,v) */
	    w = w + 1; /* addition of (u,v) increases count by 1 */
	  else if(f > 0) /* j becomes saturated by addition of (u,v) */
	    w = w + f;
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < r2) {
	  tee++;
	  f = s - aux[j];
	  if(f > 1) /* j is not saturated after addition of (u,v) */
	    w = w + 1; /* addition of (u,v) increases count by 1 */
	  else if(f > 0) /* j becomes saturated by addition of (u,v) */
	    w = w + f;
	}
      }
    }
  } else if(prop.itype == DEATH) {
    tee = aux[ix];
    if(Geyer.per) {
      /* Periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2 = dist2(u,v,x[j],y[j],Geyer.period);
	if(d2 < r2) {
	  f = s - aux[j];
	  if(f > 0) /* j is not saturated */
	    w = w + 1; /* deletion of 'ix' decreases count by 1 */
	  else {
	    f = f+1;
	    if(f > 0) {
	      /* j is not saturated after deletion of 'ix' 
		 (s must be fractional) */
	      w = w + f; 
	    }
	  }
	}
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2 = pow(u - x[j], 2) + pow(v - y[j], 2);
	if(d2 < r2) {
	  f = s - aux[j];
	  if(f > 0) /* j was not saturated */
	    w = w + 1; /* deletion of 'ix' decreases count by 1 */
	  else {
	    f = f+1; 
	    if(f > 0) {
	      /* j is not saturated after deletion of 'ix' 
		 (s must be fractional) */
	      w = w + f; 
	    }
	  }
	}
      }
    }
  } else if(prop.itype == SHIFT) { 
    /* Compute the cif at the new point, not the ratio of new/old */
    for(j=0; j<npts; j++) {
      if(j == ix) continue;
      d2 = dist2either(u,v,x[j],y[j],Geyer.period);
      if(d2 < r2) {
	tee++;
	a = aux[j];
	/* Adjust */
	dd2 = dist2either(x[ix],y[ix], x[j],y[j],Geyer.period);
	if(dd2 < r2) a = a - 1;
	b = a + 1;
	if(a < s && s < b) {
	  w = w + s - a;
	}
	else if(s >= b) w = w + 1;
      }
    }
  }

  w = w + ((tee < s) ? tee : s);

 if(Geyer.hard) {
    if(tee > 0) cifval = 0.0;
    else cifval = Geyer.beta;
  }
  else cifval = Geyer.beta*exp(Geyer.loggamma*w);
  
  return cifval;
}

void geyerupd(state, prop) 
     State state;
     Propo prop;
{
/* Declare other variables */
  int ix, npts, j, k;
  double u, v, xix, yix, r2, d2, d2old, d2new;
  double *x, *y;
  int *aux;

  aux = Geyer.aux;
  r2 = Geyer.r2;
  x = state.x;
  y = state.y;
  npts = state.npts;

  if(prop.itype == BIRTH) { 
    /* Birth */
    u = prop.u;
    v = prop.v;
    /* initialise auxiliary counter for new point */
    aux[npts] = 0; 
    /* update all auxiliary counters */
    if(Geyer.per) {
      /* periodic distance */
      for(j=0; j < npts; j++) {
	d2 = dist2(u,v,x[j],y[j],Geyer.period);
	if(d2 < r2) {
	  aux[j] += 1;
	  aux[npts] += 1;
	} 
      }
    } else {
      /* Euclidean distance */
      for(j=0; j < npts; j++) {
	d2 = pow(u-x[j], 2) + pow(v-y[j], 2);
	if(d2 < r2) {
	  aux[j] += 1;
	  aux[npts] += 1;
	} 
      }
    }
    return;
  }
  if(prop.itype == DEATH) {
    /* Death */
    ix = prop.ix;
    u = x[ix];
    v = y[ix];
    /* decrement auxiliary counter for each point */
    if(Geyer.per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	d2 = dist2(u,v,x[j],y[j],Geyer.period);
	if(d2 < r2) {
	  if(j < ix) aux[j] -= 1;
	  else aux[j-1] = aux[j] - 1;
	} else if(j >= ix) aux[j-1] = aux[j];
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	d2 = pow(u-x[j], 2) + pow(v-y[j], 2);
	if(d2 < r2) {
	  if(j < ix) aux[j] -= 1;
	  else aux[j-1] = aux[j] - 1;
	} else if(j >= ix) aux[j-1] = aux[j];
      }
    }
    return;
  }

  if(prop.itype == SHIFT) { 
    /* Shift */
    u = prop.u;
    v = prop.v;
    ix = prop.ix;
    xix = x[ix];
    yix = y[ix];
    /* recompute auxiliary counter for point 'ix' */
    aux[ix] = 0;
    /* update auxiliary counters for other points */
    if(Geyer.per) {
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2new = dist2(u,v,x[j],y[j],Geyer.period);
	d2old = dist2(xix,yix,x[j],y[j],Geyer.period);
	if(d2old >= r2 && d2new >= r2) continue;
	if(d2new < r2) {
	  /* increment neighbour count for new point */
	  aux[ix] += 1;
	  if(d2old >= r2) 
	    aux[j] += 1; /* point j gains a new neighbour */
	} else if(d2old < r2) 
	  aux[j] -= 1; /* point j loses a neighbour */
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2new = pow(u-x[j], 2) + pow(v-y[j], 2);
	d2old = pow(x[ix]-x[j], 2) + pow(y[ix]-y[j],2);
	if(d2old >= r2 && d2new >= r2) continue;
	if(d2new < r2) {
	  /* increment neighbour count for new point */
	  aux[ix] += 1;
	  if(d2old >= r2) 
	    aux[j] += 1; /* point j gains a new neighbour */
	} else if(d2old < r2) 
	  aux[j] -= 1; /* point j loses a neighbour */
      }
    }
    return;
  }
  fexitc("Unrecognised transition type; bailing out.\n");
}

Cifns GeyerCifns = { &geyerinit, &geyercif, &geyerupd, FALSE};
