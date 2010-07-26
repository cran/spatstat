#include <math.h>

/*

  Ediggatsti.c

  $Revision: 1.2 $     $Date: 2010/07/15 13:09:43 $

  C implementation of 'eval' for DiggleGatesStibbard interaction 

  Assumes point patterns are sorted in increasing order of x coordinate

*/

void Ediggatsti(nnsource, xsource, ysource, idsource, 
		nntarget, xtarget, ytarget, idtarget, 
		rrho, values) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget;
     int *idsource, *idtarget;
     double *rrho;
     /* output */
     double *values;
{
  int nsource, ntarget, j, i, ileft, iright, idsourcej;
  double xsourcej, ysourcej, xleft, xright, dx, dy, d2;
  double rho, rho2, product;

  nsource = *nnsource;
  ntarget = *nntarget;
  rho   = *rrho;

  rho2   = rho * rho;

  if(nsource == 0 || ntarget == 0) 
    return;

  ileft = iright = 0;

  for(j = 0; j < nsource; j++) {
    product = 1;
    xsourcej = xsource[j];
    ysourcej = ysource[j];
    idsourcej = idsource[j];
    /* search all points with x in [xleft, xright] */
    xleft  = xsourcej - rho;
    xright = xsourcej + rho;

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
      if(idtarget[i] != idsourcej) {
        /* squared interpoint distance */
	dx = xtarget[i] - xsourcej;
	dy = ytarget[i] - ysourcej;
	d2= dx * dx + dy * dy;
	if(d2 <= rho2) 
	  product *= sin(M_PI_2 * sqrt(d2)/rho);
      }
    }
    values[j] = log(product * product);
  }
}



