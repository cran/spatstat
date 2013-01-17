#include <R.h>
#include <math.h>
#include <stdlib.h>
#include "methas.h"
#include "dist2.h"

void fexitc(const char *msg);

/*
  Conditional intensity function for a Geyer saturation process.  
*/

typedef struct Geyer {
  /* model parameters */
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

Cdata *geyerinit(state, model, algo)
     State state;
     Model model;
     Algor algo;
{
  int i, j;
  double d2;
  Geyer *geyer;

  geyer = (Geyer *) R_alloc(1, sizeof(Geyer));

  /* Interpret model parameters*/
  geyer->gamma  = model.ipar[0];
  geyer->r      = model.ipar[1]; /* not squared any more */
  geyer->s      = model.ipar[2]; 
  geyer->r2     = geyer->r * geyer->r;
#ifdef MHDEBUG
  Rprintf("Initialising Geyer gamma=%lf, r=%lf, sat=%lf\n",
	  geyer->gamma, geyer->r, geyer->s);
#endif
  /* is the model numerically equivalent to hard core ? */
  geyer->hard   = (geyer->gamma < DOUBLE_EPS);
  geyer->loggamma = (geyer->hard) ? 0 : log(geyer->gamma);
  /* periodic boundary conditions? */
  geyer->period = model.period;
  geyer->per    = (model.period[0] > 0.0);
  /* allocate storage for auxiliary counts */
  geyer->aux = (int *) R_alloc((size_t) state.npmax, sizeof(int));
  /* Initialise auxiliary counts */
  for(i = 0; i < state.npmax; i++) 
    geyer->aux[i] = 0;

  for(i = 0; i < state.npts; i++) {
    for(j = 0; j < state.npts; j++) {
      d2 = dist2either(state.x[i], state.y[i], state.x[j], state.y[j], 
		       geyer->period);
      if(d2 < geyer->r2)
	geyer->aux[i] += 1;
    }
  }

  return((Cdata *) geyer);
}

double geyercif(prop, state, cdata)
     Propo prop;
     State state;
     Cdata *cdata;
{
  int ix, j, npts, tee;
  double u, v, d2, r2, s;
  double w, a, dd2, b, f, cifval;
  double *x, *y;
  int *aux;
  double *period;
  Geyer *geyer;
  
  geyer = (Geyer *) cdata;

  npts = state.npts;
  if(npts==0) return ((double) 1.0);

  x = state.x;
  y = state.y;
  u = prop.u;
  v = prop.v;
  ix = prop.ix;

  r2     = geyer->r2;
  s      = geyer->s;
  period = geyer->period;
  aux    = geyer->aux;

  /* 
     tee = neighbour count at the point in question;
     w   = sum of changes in (saturated) neighbour counts at other points 
  */
  tee = w = 0.0;

  if(prop.itype == BIRTH) {
    if(geyer->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	IF_CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2) {
	  tee++;
	  f = s - aux[j];
	  if(f > 1) /* j is not saturated after addition of (u,v) */
	    w = w + 1; /* addition of (u,v) increases count by 1 */
	  else if(f > 0) /* j becomes saturated by addition of (u,v) */
	    w = w + f;
	}
	END_IF_CLOSE_PERIODIC_D2
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	IF_CLOSE_D2(u,v,x[j],y[j],r2,d2) {
	  tee++;
	  f = s - aux[j];
	  if(f > 1) /* j is not saturated after addition of (u,v) */
	    w = w + 1; /* addition of (u,v) increases count by 1 */
	  else if(f > 0) /* j becomes saturated by addition of (u,v) */
	    w = w + f;
	}
	END_IF_CLOSE_D2
      }
    }
  } else if(prop.itype == DEATH) {
    tee = aux[ix];
    if(geyer->per) {
      /* Periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	IF_CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2) {
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
	END_IF_CLOSE_PERIODIC_D2
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	IF_CLOSE_D2(u,v,x[j],y[j],r2,d2) {
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
	END_IF_CLOSE_D2
      }
    }
  } else if(prop.itype == SHIFT) { 
    /* Compute the cif at the new point, not the ratio of new/old */
    if(geyer->per) {
      /* Periodic distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	IF_CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2) {
	  tee++;
	  a = aux[j];
	  /* Adjust */
	  dd2 = dist2(x[ix],y[ix], x[j],y[j],period);
	  if(dd2 < r2) a = a - 1;
	  b = a + 1;
	  if(a < s && s < b) {
	    w = w + s - a;
	  }
	  else if(s >= b) w = w + 1;
	}
	END_IF_CLOSE_PERIODIC_D2
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	IF_CLOSE_D2(u,v,x[j],y[j],r2,d2) {
	  tee++;
	  a = aux[j];
	  /* Adjust */
	  dd2 = hypot(x[ix] - x[j], y[ix] - y[j]);
	  if(dd2 < r2) a = a - 1;
	  b = a + 1;
	  if(a < s && s < b) {
	    w = w + s - a;
	  }
	  else if(s >= b) w = w + 1;
	}
	END_IF_CLOSE_PERIODIC_D2
	  }
    }
  }

  w = w + ((tee < s) ? tee : s);

 if(geyer->hard) {
    if(tee > 0) cifval = 0.0;
    else cifval = 1.0;
  }
  else cifval = exp(geyer->loggamma*w);
  
  return cifval;
}

void geyerupd(state, prop, cdata) 
     State state;
     Propo prop;
     Cdata *cdata;
{
/* Declare other variables */
  int ix, npts, j;
  double u, v, xix, yix, r2, d2, d2old, d2new;
  double *x, *y;
  int *aux;
  double *period;
  Geyer *geyer;

  geyer = (Geyer *) cdata;
  period = geyer->period;
  aux = geyer->aux;
  r2 = geyer->r2;

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
    if(geyer->per) {
      /* periodic distance */
      for(j=0; j < npts; j++) {
	IF_CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2) {
	  aux[j] += 1;
	  aux[npts] += 1;
	} 
	END_IF_CLOSE_PERIODIC_D2
      }
    } else {
      /* Euclidean distance */
      for(j=0; j < npts; j++) {
	IF_CLOSE_D2(u,v,x[j],y[j],r2,d2) {
	  aux[j] += 1;
	  aux[npts] += 1;
	} 
	END_IF_CLOSE_D2
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
    if(geyer->per) {
      /* periodic distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	IF_CLOSE_PERIODIC_D2(u,v,x[j],y[j],period,r2,d2) {
	  if(j < ix) aux[j] -= 1;
	  else aux[j-1] = aux[j] - 1;
	} else if(j >= ix) aux[j-1] = aux[j];
	END_IF_CLOSE_PERIODIC_D2
      }
    } else {
      /* Euclidean distance */
      for(j=0; j<npts; j++) {
	if(j==ix) continue;
	IF_CLOSE_D2(u,v,x[j],y[j],r2,d2) {
	  if(j < ix) aux[j] -= 1;
	  else aux[j-1] = aux[j] - 1;
	} else if(j >= ix) aux[j-1] = aux[j];
	END_IF_CLOSE_D2
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
    if(geyer->per) {
      for(j=0; j<npts; j++) {
	if(j == ix) continue;
	d2new = dist2(u,v,x[j],y[j],geyer->period);
	d2old = dist2(xix,yix,x[j],y[j],geyer->period);
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
	d2new = hypot(u-x[j],     v-y[j]);
	d2old = hypot(x[ix]-x[j], y[ix]-y[j]);
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
