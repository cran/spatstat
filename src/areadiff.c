/*

  areadiff.c

  Area difference function

  $Revision: 1.3 $ $Date: 2008/02/01 16:36:24 $

  A(u,x) = area of disc b(u,r) not covered by discs b(x_i,r) for x_i in x

  For use in area-interaction process

  Area estimated by point-counting on a fine grid

*/

#undef DEBUG

#ifdef DEBUG
#include <stdio.h>
#endif

void
areadiff(ux,uy,rad,x,y,nn,ngrid,answer) 
     double *ux, *uy;  /* coordinates of point u */
     double *rad;      /* radius */
     double *x, *y;    /* coordinate vectors for point pattern */
     int    *nn;       /* length of vectors x and y */
     int    *ngrid;    /* dimensions of point-counting grid */
     double *answer;   /* computed area */
{
  double xu, yu, dx, dy, xg, yg, r, r2, xdif, ydif;
  int i, j, k, m, jmax, n, count, found;
  xu = *ux;
  yu = *uy;
  r  = *rad;
  r2 = r * r;
  n  = *nn;
  m  = *ngrid;
  dx = dy = 2 * r / (m-1);

  /* shift point u to the origin:
     subtract (xu,yu) from all coordinate vectors */
  if(xu != 0.0 || yu != 0.0) {
#ifdef DEBUG
    fprintf(stderr, "Shifting to the origin.\n");
#endif
    for(k = 0; k < n; k++) {
      x[k] = x[k] - xu;
      y[k] = y[k] - yu;
    }
  }

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


