/*

  Efiksel.c

  $Revision: 1.1 $     $Date: 2010/07/15 04:30:39 $

  C implementation of 'eval' for Fiksel interaction (non-hardcore part)

  Assumes point patterns are sorted in increasing order of x coordinate

*/

#define OK 0
#define OVERFLOW 1

double sqrt(), exp();

void Efiksel(nnsource, xsource, ysource, 
	     nntarget, xtarget, ytarget, 
	     rrmax, kkappa, values) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget, *rrmax, *kkappa;
     /* output */
     double *values;
{
  int nsource, ntarget, j, i, ileft, iright;
  double xsourcej, ysourcej, xleft, xright, dx, dy, d2;
  double rmax, r2max, kappa, total;

  nsource = *nnsource;
  ntarget = *nntarget;
  rmax = *rrmax;
  kappa = *kkappa;

  r2max = rmax * rmax;

  if(nsource == 0 || ntarget == 0) 
    return;

  ileft = iright = 0;

  for(j = 0; j < nsource; j++) {
    total = 0;
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
	total += exp(- kappa * sqrt(d2));
    }
    values[j] = total;
  }
}



