#include <math.h>
# 1 "Kborder.cpp"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "Kborder.cpp"
# 30 "Kborder.cpp"
/* 

  Kwborder.c

  Efficient border-corrected estimate of K 
  for large datasets
  
  Kwborder()  Estimates Kinhom function.
              Requires (x,y) data to be sorted in ascending order of x
	      Expects r values to be equally spaced and starting at zero

 */


/*

  source file: 

  Kborder.cpp

  $Revision: 1.3 $     $Date: 2007/10/26 14:52:01 $

*/

double sqrt();


void Kwborder(nxy, x, y, w, b, nr, rmax, numer, denom)
     /* inputs */
     int *nxy, *nr;
     double *x, *y, *b, *rmax;
     double *w;
     /* outputs */
     double *numer, *denom;
# 72 "Kborder.cpp"
{
  int i, j, l, n, nt, n1, nt1, lmin, lmax, lup;
  double dt, tmax, tmax2, xi, yi, bi, bi2, maxsearch;
  double xleft, xright, bratio, dratio, dij, dij2, dx, dy;


  double wi, wj, wij;
# 90 "Kborder.cpp"
  n = *nxy;
  nt = *nr;

  n1 = n - 1;
  nt1 = nt - 1;

  dt = (*rmax)/(nt-1);
  tmax = *rmax;
  tmax2 = tmax * tmax;

  /* initialise */
  for(l = 0; l < nt; l++)
    numer[l] = denom[l] = 0.0;

  if(n == 0)
    return;

  for(i = 0; i < n; i++) {
    /*  --------   DENOMINATOR  -------------*/
    bi = b[i];

    wi = w[i];

    /* increment denominator for all r < b[i] */
    bratio = bi/dt;
    lmax = (int) floor(bratio);
    lup = (int) ceil(bratio);
    if(lmax == lup) --lmax;
    /* lmax is the largest integer STRICTLY less than bratio */
    lmax = (lmax <= nt1) ? lmax : nt1;
    /* increment entries 0 to lmax */
    if(lmax >= 0) {
      for(l = 0; l <= lmax; l++)
 denom[l] += wi;
    }

    /*  ----------  NUMERATOR -----------*/
    /* scan through points (x[j],y[j]) */
    bi2 = bi * bi;
    xi = x[i];
    yi = y[i];
    maxsearch = (bi < tmax) ? bi : tmax;
    xleft = xi - maxsearch;
    xright = xi + maxsearch;

    /* 
       scan backward from i-1 
       until x < xleft or until we run out 
    */
    if(i > 0) {
      for(j=i-1; x[j] >= xleft && j >= 0; j--) {
 /* squared interpoint distance */
 dx = x[j] - xi;
 dy = y[j] - yi;
 dij2= dx * dx + dy * dy;
 if(dij2 < tmax2 && dij2 < bi2) {

   wj = w[j];

   /* increment numerator for all r such that dij <= r < bi */
   dij = (double) sqrt(dij2);
   dratio = dij/dt;
   /* smallest integer greater than or equal to dratio */
   lmin = (int) ceil(dratio);
   /* increment entries lmin to lmax inclusive */
   if(lmax >= lmin) {

     wij = wi * wj;

     for(l = lmin; l <= lmax; l++)
       numer[l] += wij;
   }
 }
      }
    }

    /* 
       scan forward from i+1 
       until x > xright or until we run out 
    */
    if(i < n1) {
      for(j=i+1; x[j] <= xright && j < n; j++) {
 /* squared interpoint distance */
 dx = x[j] - xi;
 dy = y[j] - yi;
 dij2= dx * dx + dy * dy;
 if(dij2 < tmax2 && dij2 < bi2) {

   wj = w[j];

   /* increment numerator for all r such that dij <= r < bi */
   dij = (double) sqrt(dij2);
   dratio = dij/dt;
   /* smallest integer greater than or equal to dratio */
   lmin = (int) ceil(dratio);
   /* increment entries lmin to lmax inclusive */
   if(lmax >= lmin) {

     wij = wi * wj;

     for(l = lmin; l <= lmax; l++)
       numer[l] += wij;
   }
 }
      }
    }
  }
}
