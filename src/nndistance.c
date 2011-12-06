/*

  nndistance.c

  Nearest Neighbour Distances between points

  $Revision: 1.3 $     $Date: 2011/11/20 04:04:42 $

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  nndistsort    Nearest neighbour distances 
  nnwhichsort   Nearest neighbours and their distances
  nnXwhich      Nearest neighbour from one list to another
  nnXexclude    Nearest neighbour from one list to another, with overlaps

  knndsort      k-th nearest neighbour distances
  knnwhichsort  k-th nearest neighbours and their distances
*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

double sqrt();

/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

void nndistsort(n, x, y, nnd, huge)
     /* inputs */
     int *n;
     double *x, *y, *huge;
     /* output */
     double *nnd;
{ 
  int npoints, i, left, right;
  double dmin, d2, d2min, xi, yi, dx, dy, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  for(i = 0; i < npoints; i++) {
    
    R_CheckUserInterrupt();

#ifdef SPATSTAT_DEBUG
    Rprintf("\ni=%d\n", i); 
#endif

    dmin = hu;
    d2min = hu2;
    xi = x[i];
    yi = y[i];
    /* search backward */
    for(left = i - 1;
        left >= 0 && (dy = (yi - y[left])) < dmin ;
	--left)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("L");
#endif

	dx = x[left] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dy = (y[right] - yi)) < dmin ;
	++right)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("R");
#endif
	dx = x[right] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2min) {
	  d2min = d2;
	  dmin = sqrt(d2);
	}
      }
#ifdef SPATSTAT_DEBUG
    Rprintf("\n");
#endif

    nnd[i] = dmin;
  }

}

/* nnwhichsort: same as nndistsort, 
   but also returns id of nearest neighbour 
*/

void nnwhichsort(n, x, y, nnd, nnwhich, huge)
     /* inputs */
     int *n;
     double *x, *y, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints, i, left, right, which;
  double dmin, d2, d2min, xi, yi, dx, dy, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;

  for(i = 0; i < npoints; i++) {

    R_CheckUserInterrupt();

    dmin = hu;
    d2min = hu2;
    which = -1;
    xi = x[i];
    yi = y[i];
    /* search backward */
    if(i > 0){
      for(left = i - 1;
	  left >= 0 && (dy = (yi - y[left])) < dmin ;
	  --left)
	{
	  dx = x[left] - xi;
	  d2 =  dx * dx + dy * dy;
	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
	    which = left;
	  }
	}
    }

    /* search forward */
    if(i < npoints - 1) {
      for(right = i + 1;
	  right < npoints && (dy = (y[right] - yi)) < dmin ;
	  ++right)
	{
	  dx = x[right] - xi;
	  d2 =  dx * dx + dy * dy;
	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
	    which = right;
	  }
	}
    }
    nnd[i] = dmin;
    nnwhich[i] = which;
  }
}


/* 
   nnXwhich:  for TWO point patterns X and Y,
              find the nearest neighbour 
	      (from each point of X to the nearest point of Y)
	      returning both the distance and the identifier

   Requires both patterns to be sorted in order of increasing y coord
*/

void nnXwhich(n1, x1, y1, n2, x2, y2, nnd, nnwhich, huge)
     /* inputs */
     int *n1, *n2;
     double *x1, *y1, *x2, *y2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints1, npoints2, i, jleft, jright, jwhich, lastjwhich;
  double dmin, d2, d2min, x1i, y1i, dx, dy, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;

  for(i = 0; i < npoints1; i++) {

    R_CheckUserInterrupt();

    dmin = hu;
    d2min = hu2;
    jwhich = -1;
    x1i = x1[i];
    y1i = y1[i];

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {
      for(jleft = lastjwhich - 1;
	  jleft >= 0 && (dy = (y1i - y2[jleft])) < dmin ;
	  --jleft)
	{
	  dx = x2[jleft] - x1i;
	  d2 =  dx * dx + dy * dy;
	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
	    jwhich = jleft;
	  }
	}
    }

    /* search forward from previous nearest neighbour  */
    if(lastjwhich < npoints2) {
      for(jright = lastjwhich;
	  jright < npoints2 && (dy = (y2[jright] - y1i)) < dmin ;
	  ++jright)
	{
	  dx = x2[jright] - x1i;
	  d2 =  dx * dx + dy * dy;
	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
	    jwhich = jright;
	  }
	}
    }
    nnd[i] = dmin;
    nnwhich[i] = jwhich;
    lastjwhich = jwhich;
  }
}

/* 
   nnXexclude:  similar to nnXwhich
              but allows X and Y to include common points
	      (which are not to be counted as neighbours)

   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

   Requires both patterns to be sorted in order of increasing y coord
*/

void nnXexclude(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge)
     /* inputs */
     int *n1, *n2, *id1, *id2;
     double *x1, *y1, *x2, *y2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int npoints1, npoints2, i, jleft, jright, jwhich, lastjwhich, id1i;
  double dmin, d2, d2min, x1i, y1i, dx, dy, hu, hu2;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;

  for(i = 0; i < npoints1; i++) {

    R_CheckUserInterrupt();

    dmin = hu;
    d2min = hu2;
    jwhich = -1;
    x1i = x1[i];
    y1i = y1[i];
    id1i = id1[i];

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {
      for(jleft = lastjwhich - 1;
	  jleft >= 0 && (dy = (y1i - y2[jleft])) < dmin ;
	  --jleft)
	{
	  /* do not compare identical points */
	  if(id2[jleft] != id1i) {
	    dx = x2[jleft] - x1i;
	    d2 =  dx * dx + dy * dy;
	    if (d2 < d2min) {
	      d2min = d2;
	      dmin = sqrt(d2);
	      jwhich = jleft;
	    }
	  }
	}
    }

    /* search forward from previous nearest neighbour  */
    if(lastjwhich < npoints2) {
      for(jright = lastjwhich;
	  jright < npoints2 && (dy = (y2[jright] - y1i)) < dmin ;
	  ++jright)
	{
	  if(id2[jright] != id1i) {
	    dx = x2[jright] - x1i;
	    d2 =  dx * dx + dy * dy;
	    if (d2 < d2min) {
	      d2min = d2;
	      dmin = sqrt(d2);
	      jwhich = jright;
	    }
	  }
	}
    }
    nnd[i] = dmin;
    nnwhich[i] = jwhich;
    lastjwhich = jwhich;
  }
}

/* 
   knndsort 

   nearest neighbours 1:kmax

*/

void knndsort(n, kmax, x, y, nnd, huge)
     /* inputs */
     int *n, *kmax;
     double *x, *y, *huge;
     /* output matrix (npoints * kmax) in ROW MAJOR order */
     double *nnd;
{ 
  int npoints, nk, nk1, i, k, k1, left, right, unsorted;
  double d2, dminK, d2minK, xi, yi, dx, dy, hu, hu2, tmp, tmp2;
  double *dmin, *d2min;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances
     for the current point
  */

  dmin = (double *) R_alloc((size_t) nk, sizeof(double));
  d2min = (double *) R_alloc((size_t) nk, sizeof(double));

  /* loop over points */

  for(i = 0; i < npoints; i++) {

    R_CheckUserInterrupt();

#ifdef SPATSTAT_DEBUG
    Rprintf("\ni=%d\n", i); 
#endif

    /* initialise nn distances */

    dminK  = hu;
    d2minK = hu2;
    for(k = 0; k < nk; k++) {
      dmin[k] = hu;
      d2min[k] = hu2;
    }

    xi = x[i];
    yi = y[i];

    /* search backward */
    for(left = i - 1;
        left >= 0 && (dy = (yi - y[left])) < dminK ;
	--left)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("L");
#endif

	dx = x[left] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2minK) {
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  /* bubble sort */
	  unsorted = TRUE;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(dmin[k] < dmin[k1]) {
	      /* swap entries */
	      tmp = dmin[k1];
	      tmp2 = d2min[k1];
	      dmin[k1] = dmin[k];
	      d2min[k1] = d2min[k];
	      dmin[k] = tmp;
	      d2min[k] = tmp2;
	    } else {
	      unsorted = FALSE;
	    }
	  }
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dy = (y[right] - yi)) < dminK ;
	++right)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("R");
#endif
	dx = x[right] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2minK) {
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  /* bubble sort */
	  unsorted = TRUE;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(dmin[k] < dmin[k1]) {
	      /* swap entries */
	      tmp = dmin[k1];
	      tmp2 = d2min[k1];
	      dmin[k1] = dmin[k];
	      d2min[k1] = d2min[k];
	      dmin[k] = tmp;
	      d2min[k] = tmp2;
	    } else {
	      unsorted = FALSE;
	    }
	  }
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

#ifdef SPATSTAT_DEBUG
    Rprintf("\n");
#endif

    /* copy nn distances for point i 
       to output matrix in ROW MAJOR order
    */
    for(k = 0; k < nk; k++) {
      nnd[nk * i + k] = dmin[k];
    }
  }

}

/* 
   knnwhichsort 

   nearest neighbours 1:kmax

   returns distances and indices

*/

void knnwhichsort(n, kmax, x, y, nnd, nnwhich, huge)
     /* inputs */
     int *n, *kmax;
     double *x, *y, *huge;
     /* output matrices (npoints * kmax) in ROW MAJOR order */
     double *nnd;
     int    *nnwhich;
{ 
  int npoints, nk, nk1, i, k, k1, left, right, unsorted, itmp;
  double d2, dminK, d2minK, xi, yi, dx, dy, hu, hu2, tmp, tmp2;
  double *dmin, *d2min; 
  int *which;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  dmin = (double *) R_alloc((size_t) nk, sizeof(double));
  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
  which = (int *) R_alloc((size_t) nk, sizeof(int));

  /* loop over points */

  for(i = 0; i < npoints; i++) {

    R_CheckUserInterrupt();

#ifdef SPATSTAT_DEBUG
    Rprintf("\ni=%d\n", i); 
#endif

    /* initialise nn distances and indices */

    dminK  = hu;
    d2minK = hu2;
    for(k = 0; k < nk; k++) {
      dmin[k] = hu;
      d2min[k] = hu2;
      which[k] = -1;
    }

    xi = x[i];
    yi = y[i];

    /* search backward */
    for(left = i - 1;
        left >= 0 && (dy = (yi - y[left])) < dminK ;
	--left)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("L");
#endif

	dx = x[left] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2minK) {
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  which[nk1] = left;
	  /* bubble sort */
	  unsorted = TRUE;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(dmin[k] < dmin[k1]) {
	      /* swap entries */
	      tmp = dmin[k1];
	      tmp2 = d2min[k1];
	      itmp = which[k1];
	      dmin[k1] = dmin[k];
	      d2min[k1] = d2min[k];
	      which[k1] = which[k];
	      dmin[k] = tmp;
	      d2min[k] = tmp2;
	      which[k] = itmp;
	    } else {
	      unsorted = FALSE;
	    }
	  }
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dy = (y[right] - yi)) < dminK ;
	++right)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("R");
#endif
	dx = x[right] - xi;
	d2 =  dx * dx + dy * dy;
	if (d2 < d2minK) {
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  which[nk1] = right;
	  /* bubble sort */
	  unsorted = TRUE;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(dmin[k] < dmin[k1]) {
	      /* swap entries */
	      tmp = dmin[k1];
	      tmp2 = d2min[k1];
	      itmp = which[k1];
	      dmin[k1] = dmin[k];
	      d2min[k1] = d2min[k];
	      which[k1] = which[k];
	      dmin[k] = tmp;
	      d2min[k] = tmp2;
	      which[k] = itmp;
	    } else {
	      unsorted = FALSE;
	    }
	  }
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

#ifdef SPATSTAT_DEBUG
    Rprintf("\n");
#endif

    /* copy nn distances for point i 
       to output matrix in ROW MAJOR order
    */
    for(k = 0; k < nk; k++) {
      nnd[nk * i + k] = dmin[k];
      nnwhich[nk * i + k] = which[k];
    }
  }
}

