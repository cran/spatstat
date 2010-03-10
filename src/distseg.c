/*
       distseg.c

       Distance transform of a line segment pattern
       
       $Revision: 1.3 $ $Date: 2007/05/10 17:19:44 $

       Author: Adrian Baddeley

*/

#include <math.h>

void
distmap2segs(xp, yp, npoints, x0, y0, x1, y1, nsegments, epsilon, dist2, index)
     /* input */
     double	*xp, *yp;		/* pixel coordinates */
     int	*npoints;
     double	*x0, *y0, *x1, *y1;	/* line segment endpoints */
     int	*nsegments;
     double     *epsilon;               /* tolerance for short segments */
     /* output */
     double	*dist2;		        /* squared distance from pixel 
                                        to nearest line segment */
     int	*index;		        /* which line segment is closest */
{
  int	i,j, np, nseg;
  double dx,dy,leng,co,si;  /* parameters of segment */
  double xdif0,ydif0,xdif1,ydif1,xpr,ypr; /* vectors */
  double dsq0,dsq1,dsq,dsqperp; /* squared distances */
  double eps;

  np   = *npoints;
  nseg = *nsegments;
  eps  = *epsilon;
  for(j = 0; j < nseg; j++) {
    dx = x1[j] - x0[j];
    dy = y1[j] - y0[j];
    leng = sqrt(dx*dx + dy*dy);
    co = dx/leng;
    si = dy/leng;
    for(i = 0; i < np; i++) {
      /* vectors from pixel to segment endpoints */
      xdif0 =  xp[i] - x0[j];
      ydif0 =  yp[i] - y0[j];
      xdif1 =  xp[i] - x1[j];
      ydif1 =  yp[i] - y1[j];
      /* squared distances to segment endpoints */
      dsq0 = xdif0*xdif0 + ydif0*ydif0;
      dsq1 = xdif1*xdif1 + ydif1*ydif1;
      dsq = (dsq0 < dsq1) ? dsq0 : dsq1;
      if(leng > eps) {
	/* rotate pixel around 1st endpoint of segment
	   so that line segment lies in x axis */
	xpr = xdif0 * co + ydif0 * si;
	ypr = -xdif0 * si + ydif0 * co;
	/* perpendicular distance applies only in perpendicular region */
	if(xpr >= 0.0 && xpr <= leng) {
	  dsqperp = ypr*ypr;
	  if(dsqperp < dsq) dsq = dsqperp;
	}
      }
      if(j == 0) {
	dist2[i] = dsq;
	index[i] = 0;
      } else if(dist2[i] > dsq) {
	dist2[i] = dsq;
	index[i] = j;
      }
    }
  }
}	
