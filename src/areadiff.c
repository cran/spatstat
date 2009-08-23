/*

  areadiff.c

  Area difference function

  $Revision: 1.7 $ $Date: 2009/08/12 00:20:28 $

  A(x,r) = area of disc b(0,r) not covered by discs b(x_i,r) for x_i in x
  
  Area estimated by point-counting on a fine grid

  For use in area-interaction model and related calculations

*/

#undef DEBUG

#ifdef DEBUG
#include <stdio.h>
#endif

/* 
   Original version areadiff()

   1 point u

   No trimming of discs

*/

void
areadiff(rad,x,y,nn,ngrid,answer) 
     /* inputs */
     double *rad;      /* radius */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nn;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     /* output */
     double *answer;   /* computed area */
{
  double dx, dy, xg, yg, r, r2, xdif, ydif;
  int i, j, k, m, n, count, found;
  r  = *rad;
  r2 = r * r;
  n  = *nn;
  m  = *ngrid;
  dx = dy = 2 * r / (m-1);

  count = 0;

  /* run through grid points */
  for(i = 0, xg = -r; i < m; i++, xg += dx) 
    for(j = 0, yg = -r; j < m; j++, yg += dy)
      /* test for inside disc */
      if(xg * xg + yg * yg < r2) {
#ifdef DEBUG
	fprintf(stderr, "\n\n (xg,yg) = (%lf, %lf)\n", xg, yg);
#endif
	/* run through data points seeking one close to (xy, yg) */
	found = 0;
	if(n > 0) 
	  for(k = 0; k < n && found == 0; k++) {
#ifdef DEBUG
	    fprintf(stderr, "(x[%d],y[%d]) = (%lf,%lf)\n", k, k, x[k], y[k]);
#endif
	    xdif = x[k] - xg;
	    ydif = y[k] - yg;
	    if(xdif * xdif + ydif * ydif < r2) {
	      found = 1;
#ifdef DEBUG
	      fprintf(stderr, "(x[%d], y[%d]) = (%lf, %lf) covers!\n", 
		      k, k, x[k], y[k]);
#endif
	    }
	  }
	if(found == 0) {
	  ++count;
#ifdef DEBUG
	  fprintf(stderr, "---------------Incrementing count\n");
#endif
	    }
      }

#ifdef DEBUG
  fprintf(stderr, "Count = %d\n", count);
#endif
  
  /* calculate area */
  *answer = ((double) count) * dx * dy;
}


/* similar function, handles multiple values of 'r' */

void
areadifs(rad,nrads,x,y,nxy,ngrid,answer) 
     /* inputs */
     double *rad;      /* vector of radii */
     int    *nrads;     /* length of 'rads' */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nxy;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     /* output */
     double *answer;   /* computed areas (vector of length 'nrads') */
{
  double dx, dy, xg, yg, r, r2, xdif, ydif;
  int i, j, k, l, m, n, nr, count, found;

  n  = *nxy;
  nr = *nrads;
  m  = *ngrid;

  /* run through radii */
  for(l = 0; l < nr; l++) {
    r  = rad[l];
    if(r == 0.0) {
      answer[l] = 0.0;
    } else {
      r2 = r * r;
      dx = dy = 2 * r / (m-1);
      count = 0;

      /* run through grid points */
      for(i = 0, xg = -r; i < m; i++, xg += dx) 
	for(j = 0, yg = -r; j < m; j++, yg += dy)
	  /* test for inside disc */
	  if(xg * xg + yg * yg < r2) {
#ifdef DEBUG
	    fprintf(stderr, "\n\n (xg,yg) = (%lf, %lf)\n", xg, yg);
#endif
	    /* run through data points seeking one close to (xy, yg) */
	    found = 0;
	    if(n > 0)
	      for(k = 0; k < n && found == 0; k++) {
#ifdef DEBUG
		fprintf(stderr, "(x[%d],y[%d]) = (%lf,%lf)\n", k, k, x[k], y[k]);
#endif
		xdif = x[k] - xg;
		ydif = y[k] - yg;
		if(xdif * xdif + ydif * ydif < r2) {
		  found = 1;
#ifdef DEBUG
		  fprintf(stderr, "(x[%d], y[%d]) = (%lf, %lf) covers!\n", 
			  k, k, x[k], y[k]);
#endif
		}
	      }
	    if(found == 0) {
	      ++count;
#ifdef DEBUG
	      fprintf(stderr, "---------------Incrementing count\n");
#endif
	    }
	  }
      /* end of loop through data points */

#ifdef DEBUG
      fprintf(stderr, "Count = %d\n", count);
#endif
  
      /* calculate area for this value of r*/
      answer[l] = ((double) count) * dx * dy;
    }
    /* end of if(r==0).. else {...} */
  }
  /* end of loop over r values */
}


/*
    Modified version

    multiple test points u
    
    discs constrained inside a rectangle

*/

void
areaBdif(rad,nrads,x,y,nxy,ngrid,x0,y0,x1,y1,answer) 
     /* inputs */
     double *rad;      /* vector of radii */
     int    *nrads;     /* length of 'rads' */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nxy;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     double *x0,*y0,*x1,*y1;  /* constraint rectangle */
     /* output */
     double *answer;   /* computed areas (vector of length 'nrads') */
{
  double dx, dy, xg, yg, r, r2, xdif, ydif, xmin, ymin, xmax, ymax;
  int i, j, k, l, m, n, nr, count, found;

  n  = *nxy;
  nr = *nrads;
  m  = *ngrid;

  xmin = *x0;
  ymin = *y0;
  xmax = *x1;
  ymax = *y1;

  /* run through radii */
  for(l = 0; l < nr; l++) {
    r  = rad[l];
    if(r == 0.0) {
      answer[l] = 0.0;
    } else {
      r2 = r * r;
      dx = dy = 2 * r / (m-1);
      count = 0;

      /* run through grid points */
      for(i = 0, xg = -r; i < m; i++, xg += dx) {
	if(xg >= xmin && xg <= xmax) {
	  for(j = 0, yg = -r; j < m; j++, yg += dy) {
	    if(yg >= ymin && yg <= ymax) {
	      /* test for inside disc */
	      if(xg * xg + yg * yg < r2) {
#ifdef DEBUG
		fprintf(stderr, "\n\n (xg,yg) = (%lf, %lf)\n", xg, yg);
#endif
	    /* run through data points seeking one close to (xy, yg) */
		found = 0;
		if(n > 0)
		  for(k = 0; k < n && found == 0; k++) {
#ifdef DEBUG
		    fprintf(stderr, "(x[%d],y[%d]) = (%lf,%lf)\n", 
			    k, k, x[k], y[k]);
#endif
		    xdif = x[k] - xg;
		    ydif = y[k] - yg;
		    if(xdif * xdif + ydif * ydif < r2) {
		      found = 1;
#ifdef DEBUG
		      fprintf(stderr, "(x[%d], y[%d]) = (%lf, %lf) covers!\n", 
			      k, k, x[k], y[k]);
#endif
		    }
		  }
		if(found == 0) {
		  ++count;
#ifdef DEBUG
		  fprintf(stderr, "---------------Incrementing count\n");
#endif
		}
	      }
	    }
	  }
	}
      }
      /* end of loop through data points */

#ifdef DEBUG
      fprintf(stderr, "Count = %d\n", count);
#endif
  
      /* calculate area for this value of r*/
      answer[l] = ((double) count) * dx * dy;
    }
    /* end of if(r==0).. else {...} */
  }
  /* end of loop over r values */
}


