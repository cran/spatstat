/*

  Estrauss.c

  $Revision: 1.2 $     $Date: 2010/07/12 13:25:02 $

  C implementation of 'eval' for Strauss interaction

  Calculates number of data points within distance r of each quadrature point
  (when 'source' = quadrature points, 'target' = data points)

  Assumes point patterns are sorted in increasing order of x coordinate

*/

#define OK 0
#define OVERFLOW 1

double sqrt();

void closepaircounts(nnsource, xsource, ysource, 
		     nntarget, xtarget, ytarget, 
		     rrmax, counts) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget, *rrmax;
     /* output */
     int *counts;
{
  int nsource, ntarget, j, i, ileft, iright, counted;
  double xsourcej, ysourcej, rmax, r2max, xleft, xright, dx, dy, d2;

  nsource = *nnsource;
  ntarget = *nntarget;
  rmax = *rrmax;
  r2max = rmax * rmax;

  if(nsource == 0 || ntarget == 0) 
    return;

  ileft = iright = 0;

  for(j = 0; j < nsource; j++) {
    counted = 0;
    xsourcej = xsource[j];
    ysourcej = ysource[j];
    /* search all points with x in [xleft, xright] */
    xleft  = xsourcej - rmax;
    xright = xsourcej + rmax;

    /* 
       adjust scope of search [ileft, iright]

    */
    while((ileft+1 < ntarget) && xtarget[ileft] < xleft)
      ++ileft;

    while((iright+1 < ntarget) && xtarget[iright+1] <= xright)
      ++iright;

    /* 
       process from ileft to iright
    */
    for(i=ileft; i <= iright; i++) {
      /* squared interpoint distance */
      dx = xtarget[i] - xsourcej;
      dy = ytarget[i] - ysourcej;
      d2= dx * dx + dy * dy;
      if(d2 <= r2max)
	++counted;
    }
    counts[j] = counted;
  }
}



