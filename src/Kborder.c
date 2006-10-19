#include <math.h>
# 1 "Kborder.cpp"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "Kborder.cpp"
# 18 "Kborder.cpp"
/*

  Kborder.c

  Efficient border-corrected estimate of K 
  for large datasets
  
  Kborder()   Estimates K function.
              Requires (x,y) data to be sorted in ascending order of x
	      Expects r values to be equally spaced and starting at zero
*/
# 44 "Kborder.cpp"
/*

  source file: 

  Kborder.cpp

  $Revision: 1.1 $     $Date: 2006/10/19 10:22:21 $

*/

double sqrt();
# 65 "Kborder.cpp"
void Kborder(nxy, x, y, b, nr, rmax, numer, denom)
     /* inputs */
     int *nxy, *nr;
     double *x, *y, *b, *rmax;
     /* outputs */
     int *numer, *denom;

{
  int i, j, k, l, n, nt, n1, nt1, lmin, lmax, lup;
  double dt, tmax, tmax2, xi, yi, bi, bi2, maxsearch;
  double xleft, xright, bratio, dratio, dij, dij2, dd, dx, dy;
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
    numer[l] = denom[l] = 0;

  if(n == 0)
    return;

  for(i = 0; i < n; i++) {
    /*  --------   DENOMINATOR  -------------*/
    bi = b[i];



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
 denom[l] += 1;
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



   /* increment numerator for all r such that dij <= r < bi */
   dij = (double) sqrt(dij2);
   dratio = dij/dt;
   /* smallest integer greater than or equal to dratio */
   lmin = (int) ceil(dratio);
   /* increment entries lmin to lmax inclusive */
   if(lmax >= lmin) {



     for(l = lmin; l <= lmax; l++)
       numer[l] += 1;
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



   /* increment numerator for all r such that dij <= r < bi */
   dij = (double) sqrt(dij2);
   dratio = dij/dt;
   /* smallest integer greater than or equal to dratio */
   lmin = (int) ceil(dratio);
   /* increment entries lmin to lmax inclusive */
   if(lmax >= lmin) {



     for(l = lmin; l <= lmax; l++)
       numer[l] += 1;
   }
 }
      }
    }
  }
}
