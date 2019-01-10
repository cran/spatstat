#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"
#include "looptest.h"
/*

  Estrauss.c

  $Revision: 1.5 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  C implementation of 'eval' for Strauss interaction

  Calculates number of data points within distance r of each quadrature point
  (when 'source' = quadrature points, 'target' = data points)

  Assumes point patterns are sorted in increasing order of x coordinate

*/

double sqrt();

void Ccrosspaircounts(nnsource, xsource, ysource, 
		     nntarget, xtarget, ytarget, 
		     rrmax, counts) 
/* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget, *rrmax;
     /* output */
     int *counts;
{
  int nsource, ntarget, maxchunk, j, i, ileft, counted;
  double xsourcej, ysourcej, rmax, r2max, r2maxpluseps, xleft, dx, dy, dx2, d2;

  nsource = *nnsource;
  ntarget = *nntarget;
  rmax = *rrmax;
  r2max = rmax * rmax;
  r2maxpluseps = r2max + EPSILON(r2max);

  if(nsource == 0 || ntarget == 0) 
    return;

  ileft = 0;

  OUTERCHUNKLOOP(j, nsource, maxchunk, 65536) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nsource, maxchunk, 65536) {
      counted = 0;
      xsourcej = xsource[j];
      ysourcej = ysource[j];
      /* 
	 adjust starting point
      */
      xleft  = xsourcej - rmax;
      while((xtarget[ileft] < xleft) && (ileft+1 < ntarget))
	++ileft;

      /* 
	 process from ileft to iright
      */
      for(i=ileft; i < ntarget; i++) {
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > r2maxpluseps)
	  break;
	dy = ytarget[i] - ysourcej;
	d2 = dx2 + dy * dy;
	if(d2 <= r2max)
	  ++counted;
      }
      counts[j] = counted;
    }
  }
}

