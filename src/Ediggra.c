/*

  Ediggra.c

  $Revision: 1.2 $     $Date: 2010/07/15 13:09:43 $

  C implementation of 'eval' for DiggleGratton interaction (exponentiated)

  Assumes point patterns are sorted in increasing order of x coordinate

*/

#define OK 0
#define OVERFLOW 1

double sqrt();

void Ediggra(nnsource, xsource, ysource, idsource, 
	     nntarget, xtarget, ytarget, idtarget, 
	     ddelta, rrho, values) 
     /* inputs */
     int *nnsource, *nntarget;
     double *xsource, *ysource, *xtarget, *ytarget;
     int *idsource, *idtarget;
     double *ddelta, *rrho;
     /* output */
     double *values;
{
  int nsource, ntarget, j, i, ileft, iright, idsourcej;
  double xsourcej, ysourcej, xleft, xright, dx, dy, d2;
  double delta, rho, delta2, rho2, rhominusdelta;
  double product;

  nsource = *nnsource;
  ntarget = *nntarget;
  delta = *ddelta;
  rho   = *rrho;

  rho2   = rho * rho;
  delta2 = delta * delta;
  rhominusdelta = rho - delta;

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
	if(d2 <= rho2) {
	  if(d2 <= delta2)
	    product = 0;
	  else 
	    product *= (sqrt(d2) - delta)/rhominusdelta;
	}
      }
    }
    /* allow for numerical errors */
    values[j] = product;
  }
}



