/*

  Egeyer.c

  $Revision: 1.2 $     $Date: 2010/07/13 05:45:46 $

  Part of C implementation of 'eval' for Geyer interaction

  Calculates change in saturated count 

  (xquad, yquad): quadscheme 
  (xdata, ydata): data
  tdata: unsaturated pair counts for data pattern
  quadtodata[j] = i   if quad[j] == data[i]  (indices start from ZERO)
  
  Assumes point patterns are sorted in increasing order of x coordinate

*/

#define OK 0
#define OVERFLOW 1

double sqrt();

void Egeyer(nnquad, xquad, yquad, quadtodata,
		 nndata, xdata, ydata, tdata,
		 rrmax, ssat, result) 
     /* inputs */
     int *nnquad, *nndata, *quadtodata, *tdata;
     double *xquad, *yquad, *xdata, *ydata, *rrmax, *ssat;
     /* output */
     double *result;
{
  int nquad, ndata, j, i, ileft, iright, total, dataindex, isdata;
  double xquadj, yquadj, rmax, sat, r2max, xleft, xright, dx, dy, d2;
  double tbefore, tafter, satbefore, satafter, delta;

  nquad = *nnquad;
  ndata = *nndata;
  rmax  = *rrmax;
  sat   = *ssat;

  r2max = rmax * rmax;

  if(nquad == 0 || ndata == 0) 
    return;

  ileft = iright = 0;

  for(j = 0; j < nquad; j++) {
    total = 0;
    xquadj = xquad[j];
    yquadj = yquad[j];
    dataindex = quadtodata[j];
    isdata = (dataindex >= 0);
    /* search all points with x in [xleft, xright] */
    xleft  = xquadj - rmax;
    xright = xquadj + rmax;

    /* 
       adjust scope of search [ileft, iright]

    */
    while((ileft+1 < ndata) && xdata[ileft] < xleft)
      ++ileft;

    while((iright+1 < ndata) && xdata[iright+1] <= xright)
      ++iright;

    /* 
       process from ileft to iright
    */
    for(i=ileft; i <= iright; i++) {
      if(i != dataindex) {
	/* squared interpoint distance */
	dx = xdata[i] - xquadj;
	dy = ydata[i] - yquadj;
	d2= dx * dx + dy * dy;
	if(d2 <= r2max) {
	  /* effect of adding dummy point j or 
	     negative effect of removing data point */
	  tbefore = tdata[i];
	  tafter  = tbefore + ((isdata) ? -1 : 1);
	  /* effect on saturated values */
	  satbefore = (double) ((tbefore < sat)? tbefore : sat);
	  satafter  = (double) ((tafter  < sat)? tafter  : sat);
	  /* sum changes over all i */
	  delta = satafter - satbefore; 
	  total += ((isdata) ? -delta : delta);
	}
      }
    }
    result[j] = total;
  }
}



