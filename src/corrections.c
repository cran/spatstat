/*

  corrections.c

  Edge corrections

  $Revision: 1.5 $     $Date: 2007/10/26 14:57:53 $

 */

#undef DEBUG

#include <math.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define PI 3.1415926535898
#define YES (0 == 0)
#define NO  (0 == 1)

#define BETWEEN(X,X0,X1) (((X) - (X0)) * ((X) - (X1)) <= 0)

#define UNDER(X,Y,X0,Y0,X1,Y1) \
  (((Y1) - (Y0)) * ((X) - (X0)) >= ((Y) - (Y0)) * ((X1)- (X0)))

#define UNDERNEATH(X,Y,X0,Y0,X1,Y1) \
    (((X0) < (X1)) ? UNDER(X,Y,X0,Y0,X1,Y1) : UNDER(X,Y,X1,Y1,X0,Y0))

#define TESTINSIDE(X,Y,X0,Y0,X1,Y1) \
  (BETWEEN(X,X0,X1) && UNDERNEATH(X, Y, X0, Y0, X1, Y1))


void ripleybox(nx, x, y, rmat, nr, xmin, ymin, xmax, ymax,  epsilon, out)
     int *nx, *nr;  /* dimensions */
     double *x, *y; /* coordinate vectors of length nx */
     double *rmat;  /* matrix nx by nr  */
     double *xmin, *ymin, *xmax, *ymax;  /* box dimensions */
     double *epsilon; 
     double *out;  /* output matrix nx by nr */
{
  int i, j, n, m, ijpos, ncor;
  double xx, yy, x0, y0, x1, y1, dL, dR, dU, dD, aL, aU, aD, aR, rij;
  double cL, cU, cD, cR, bLU, bLD, bRU, bRD, bUL, bUR, bDL, bDR;
  double corner, extang;
  double eps;

  n  = *nx;
  m  = *nr;
  x0 = *xmin;
  y0 = *ymin;
  x1 = *xmax;
  y1 = *ymax;
  eps = *epsilon;
  for(i = 0; i < n; i++) {
     xx = x[i];
     yy = y[i];
  /* 
  perpendicular distance from point to each edge of rectangle
  L = left, R = right, D = down, U = up
  */
     dL = xx - x0;
     dR = x1 - xx;
     dD = yy - y0;
     dU = y1 - yy;

    /*
      test for corner of the rectangle
    */
#define ABS(X) (((X) >= 0) ? (X) : (-X))
#define SMALL(X) ((ABS(X) < eps) ? 1 : 0)

     ncor = SMALL(dL) + SMALL(dR) + SMALL(dD) + SMALL(dU);
     corner = (ncor >= 2) ? YES : NO;
  
    /* 
      angle between 
            - perpendicular to edge of rectangle
      and 
            - line from point to corner of rectangle

    */
     bLU = atan2(dU, dL);
     bLD = atan2(dD, dL);
     bRU = atan2(dU, dR);
     bRD = atan2(dD, dR);
     bUL = atan2(dL, dU);
     bUR = atan2(dR, dU);
     bDL = atan2(dL, dD);
     bDR = atan2(dR, dD);

     for(j = 0; j < m; j++) {
       ijpos = j * n + i;
       rij = rmat[ijpos];
#ifdef DEBUG
       fprintf(stderr, "rij = %lf\n", rij);
#endif
       /*
	 half the angle subtended by the intersection between
         the circle of radius r[i,j] centred on point i
         and each edge of the rectangle (prolonged to an infinite line)
       */
       aL = (dL < rij) ? acos(dL/rij) : 0.0;
       aR = (dR < rij) ? acos(dR/rij) : 0.0;
       aD = (dD < rij) ? acos(dD/rij) : 0.0;
       aU = (dU < rij) ? acos(dU/rij) : 0.0;
#ifdef DEBUG
       fprintf(stderr, "aL = %lf\n", aL);
       fprintf(stderr, "aR = %lf\n", aR);
       fprintf(stderr, "aD = %lf\n", aD);
       fprintf(stderr, "aU = %lf\n", aU);
#endif
       /* apply maxima */

       cL = MIN(aL, bLU) + MIN(aL, bLD);
       cR = MIN(aR, bRU) + MIN(aR, bRD);
       cU = MIN(aU, bUL) + MIN(aU, bUR);
       cD = MIN(aD, bDL) + MIN(aD, bDR);
#ifdef DEBUG
       fprintf(stderr, "cL = %lf\n", cL);
       fprintf(stderr, "cR = %lf\n", cR);
       fprintf(stderr, "cD = %lf\n", cD);
       fprintf(stderr, "cU = %lf\n", cU);
#endif

       /* total exterior angle over 2 pi */
       extang = (cL + cR + cU + cD)/(2 * PI);

       /* add pi/2 for corners */
       if(corner) 
	 extang += 1/4;

#ifdef DEBUG
       fprintf(stderr, "extang = %lf\n", extang);
#endif
       /* OK, now compute weight */
       out[ijpos] = 1 / (1 - extang);
     }
  }
}


void ripleypoly(nc, xc, yc, nr, rmat, nseg, x0, y0, x1, y1, out) 
     int *nc, *nr, *nseg;
     double *xc, *yc, *rmat;
     double *x0, *y0, *x1, *y1;
     double *out;
{
  int n, m, i, j, k, l, nradperpt, ncut, nchanges;
  double xcentre, ycentre, xx0, yy0, xx1, yy1, xx01, yy01;
  double x, y, radius, radius2, dx0, dx1, dy0;
  double a, b, c, t, det, sqrtdet, tmp;
  double theta[6], delta[7], tmid[7];
  double xtest, ytest, contrib, total;

  n = *nc;
  nradperpt = *nr;
  m = *nseg;

  for(i = 0; i < n; i++) {
    xcentre = xc[i];
    ycentre = yc[i];
#ifdef DEBUG
    fprintf(stderr, "centre = (%lf, %lf)\n", xcentre, ycentre);
#endif

    for(j = 0; j < nradperpt; j++) {
      radius = rmat[ j * n + i];
      radius2 = radius * radius;
#ifdef DEBUG
       fprintf(stderr, "radius = %lf\n", radius);
#endif

      total = 0.0;
      for(k=0; k < m; k++) {
#ifdef DEBUG
       fprintf(stderr, "k = %d\n", k);
#endif
	ncut = 0;
	xx0 = x0[k];
	yy0 = y0[k];
	xx1 = x1[k];
	yy1 = y1[k];
#ifdef DEBUG
       fprintf(stderr, "(%lf,%lf) to (%lf,%lf)\n", xx0, yy0, xx1, yy1);
#endif
	/* intersection with left edge */
	dx0 = xx0 - xcentre;
	det = radius2 - dx0 * dx0;
	if(det > 0) {
	  sqrtdet = sqrt(det);
	  y = ycentre + sqrtdet;
	  if(y < yy0) {
	    theta[ncut] = atan2(y - ycentre, dx0);
#ifdef DEBUG
	    fprintf(stderr, "cut left at theta= %lf\n", theta[ncut]);
#endif
	    ncut++;
	  }
	  y = ycentre - sqrtdet;
	  if(y < yy0) {
	    theta[ncut] = atan2(y-ycentre, dx0);
#ifdef DEBUG
	    fprintf(stderr, "cut left at theta= %lf\n", theta[ncut]);
#endif
	    ncut++;
	  }
	} else if(det == 0) {
	  if(ycentre < yy0) {
	    theta[ncut] = atan2(0.0, dx0);
#ifdef DEBUG
	    fprintf(stderr, "tangent left at theta= %lf\n", theta[ncut]);
#endif
	    ncut++;
	  }
	}
	/* intersection with right edge */
	dx1 = xx1 - xcentre;
	det = radius2 - dx1 * dx1;
	if(det > 0) {
	  sqrtdet = sqrt(det);
	  y = ycentre + sqrtdet;
	  if(y < yy1) {
	    theta[ncut] = atan2(y - ycentre, dx1);
#ifdef DEBUG
	    fprintf(stderr, "cut right at theta= %lf\n", theta[ncut]);
#endif
	    ncut++;
	  }
	  y = ycentre - sqrtdet;
	  if(y < yy1) {
	    theta[ncut] = atan2(y - ycentre, dx1);
#ifdef DEBUG
	    fprintf(stderr, "cut right at theta= %lf\n", theta[ncut]);
#endif
	    ncut++;
	  }
	} else if(det == 0) {
	  if(ycentre < yy1) {
	    theta[ncut] = atan2(0.0, dx1);
#ifdef DEBUG
	    fprintf(stderr, "tangent right at theta= %lf\n", theta[ncut]);
#endif
	    ncut++;
	  }
	}
	/* intersection with top segment */
	xx01 = xx1 - xx0;
	yy01 = yy1 - yy0;
	dy0  = yy0 - ycentre;
	a = xx01 * xx01 + yy01 * yy01;
	b = 2 * (xx01 * dx0 + yy01 * dy0);
	c = dx0 * dx0 + dy0 * dy0 - radius2;
	det = b * b - 4 * a * c;
	if(det > 0) {
	  sqrtdet = sqrt(det);
	  t = (sqrtdet - b)/(2 * a);
	  if(t >= 0 && t <= 1) {
	    x = xx0 + t * xx01;
	    y = yy0 + t * yy01;
	    theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUG
	    fprintf(stderr, "hits segment: t = %lf, theta = %lf\n", 
		    t, theta[ncut]);
#endif
	    ++ncut;
	  }
	  t = (-sqrtdet - b)/(2 * a);
	  if(t >= 0 && t <= 1) {
	    x = xx0 + t * xx01;
	    y = yy0 + t * yy01;
	    theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUG
	    fprintf(stderr, "hits segment: t = %lf, theta = %lf\n", 
		    t, theta[ncut]);
#endif
	    ++ncut;
	  }
	} else if(det == 0) {
	  t = - b/(2 * a);
	  if(t >= 0 && t <= 1) {
	    x = xx0 + t * xx01;
	    y = yy0 + t * yy01;
	    theta[ncut] = atan2(y - ycentre, x - xcentre);
#ifdef DEBUG
	    fprintf(stderr, "tangent to segment: t = %lf, theta = %lf\n", 
		    t, theta[ncut]);
#endif
	    ++ncut;
	  }
	}
	/* for safety, force all angles to be in range [0, 2 * pi] */
	if(ncut > 0) 
	  for(l = 0; l < ncut; l++)
	    if(theta[l] < 0) 
	      theta[l] += 2 * PI;

	/* sort angles */
	if(ncut > 1) {
	  do {
	    nchanges = 0;
	    for(l = 0; l < ncut - 1; l++) {
	      if(theta[l] > theta[l+1]) {
		/* swap */
		++nchanges;
		tmp = theta[l];
		theta[l] = theta[l+1];
		theta[l+1] = tmp;
	      }
	    }
	  } while(nchanges > 0);
	}
#ifdef DEBUG
	if(ncut > 0) {
	  for(l = 0; l < ncut; l++)
	    fprintf(stderr, "theta[%d] = %lf\n", l, theta[l]);
	}
#endif
	/* compute length of circumference inside polygon */
	if(ncut == 0) {
	  /* entire circle is either in or out */
	  xtest = xcentre + radius;
	  ytest = ycentre;
	  if(TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) 
	    contrib = 2 * PI;
	  else 
	    contrib = 0.0;
	} else {
	  /* find midpoints and lengths of pieces (adding theta = ) */
	  delta[0] = theta[0];
	  tmid[0] = theta[0]/2;
	  if(ncut > 1) {
	    for(l = 1; l < ncut; l++) {
	      delta[l] = theta[l] - theta[l-1];
	      tmid[l] = (theta[l] + theta[l-1])/2;
	    }
	  }
	  delta[ncut] = 2 * PI - theta[ncut - 1];
	  tmid[ncut] = (2 * PI + theta[ncut-1])/2;
	  contrib = 0.0;
	  for(l = 0; l <= ncut; l++) {
#ifdef DEBUG
	    fprintf(stderr, "delta[%d] = %lf\n", l, delta[l]);
#endif
	    xtest = xcentre + radius * cos(tmid[l]);
	    ytest = ycentre + radius * sin(tmid[l]);
	    if(TESTINSIDE(xtest, ytest, xx0, yy0, xx1, yy1)) {
	      contrib += delta[l];
#ifdef DEBUG 
	      fprintf(stderr, "... inside\n");
	    } else {
	      fprintf(stderr, "... outside\n");
#endif
	    }

	  }
	}
	/* multiply by sign of trapezium */
	if(xx0  < xx1)
	  contrib *= -1;

#ifdef DEBUG
	fprintf(stderr, "contrib = %lf\n", contrib);
#endif
	total += contrib;
      }
      out[ j * n + i] = total;
#ifdef DEBUG
	fprintf(stderr, "total = %lf\n", total);
#endif
    }
  }
}




