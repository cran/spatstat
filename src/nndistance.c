/*

  nndistance.c

  Nearest Neighbour Distances between points

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2012
  Licence: GNU Public Licence >= 2

  $Revision: 1.14 $     $Date: 2012/03/14 03:35:16 $

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  nndistsort    Nearest neighbour distances 
  nnwhichsort   Nearest neighbours
  nnsort        Nearest neighbours & distances

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

/* ------------------- one point pattern X --------------------- */

/* 
   nndistsort: nearest neighbour distances 
*/

#undef FNAME
#undef DIST
#undef WHICH
#define FNAME nndistsort
#define DIST
#include "nndist.h"

/* 
   nnwhichsort: id of nearest neighbour 
*/

#undef FNAME
#undef DIST
#undef WHICH
#define FNAME nnwhichsort
#define WHICH
#include "nndist.h"

/* 
   nnsort: distance & id of nearest neighbour 
*/

#undef FNAME
#undef DIST
#undef WHICH
#define FNAME nnsort
#define DIST
#define WHICH
#include "nndist.h"

/* --------------- two distinct point patterns X and Y  ----------------- */

/* 
   nnXdist:  nearest neighbour distance
	      (from each point of X to the nearest point of Y)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXdist
#define DIST
#include "nndistX.h"

/* 
   nnXwhich:  nearest neighbour id
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXwhich
#define WHICH
#include "nndistX.h"

/* 
   nnX:  nearest neighbour distance and id
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnX
#define DIST
#define WHICH
#include "nndistX.h"

/* --------------- two point patterns X and Y with common points --------- */

/*
   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].
*/

/* 
   nnXEdist:  similar to nnXdist
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXEdist
#define DIST
#define EXCLUDE
#include "nndistX.h"

/* 
   nnXEwhich:  similar to nnXwhich
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXEwhich
#define WHICH
#define EXCLUDE
#include "nndistX.h"

/* 
   nnXE:  similar to nnX
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#define FNAME nnXE
#define DIST
#define WHICH
#define EXCLUDE
#include "nndistX.h"


/*    -------------- k-th nearest neighbours ------------------ */

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
  int npoints, maxchunk, nk, nk1, i, k, k1, left, right, unsorted;
  double d2, d2minK, xi, yi, dx, dy, dy2, hu, hu2, tmp;
  double *d2min;

  hu = *huge;
  hu2 = hu * hu;

  npoints = *n;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < npoints) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints) maxchunk = npoints;

    for(; i < maxchunk; i++) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

    /* initialise nn distances */

      d2minK = hu2;
      for(k = 0; k < nk; k++)
	d2min[k] = hu2;

      xi = x[i];
      yi = y[i];

      /* search backward */
      for(left = i - 1; left >= 0; --left)
	{

#ifdef SPATSTAT_DEBUG
	Rprintf("L");
#endif
	  dy = yi - y[left];
	  dy2 = dy * dy; 
	  if(dy2 > d2minK)
	    break;

	  dx = x[left] - xi;
	  d2 =  dx * dx + dy2;
	  if (d2 < d2minK) {
	    /* overwrite last entry */
	    d2min[nk1] = d2;
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
	      } else {
		unsorted = FALSE;
	      }
	    }
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}
      
    /* search forward */
      for(right = i + 1; right < npoints; ++right)
	{

#ifdef SPATSTAT_DEBUG
	  Rprintf("R");
#endif
	  dy = y[right] - yi;
	  dy2 = dy * dy;
	  if(dy2 > d2minK) 
	    break;

	  dx = x[right] - xi;
	  d2 =  dx * dx + dy2;
	  if (d2 < d2minK) {
	    /* overwrite last entry */
	    d2min[nk1] = d2;
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp = d2min[k1];
		d2min[k1] = d2min[k];
		d2min[k] = tmp;
	      } else {
		unsorted = FALSE;
	      }
	    }
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}

      /* search finished for point i */

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif
      
      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) 
	nnd[nk * i + k] = sqrt(d2min[k]);
      
      /* end i loop */
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
  int npoints, maxchunk, nk, nk1, i, k, k1, left, right, unsorted, itmp;
  double d2, d2minK, xi, yi, dx, dy, dy2, hu, hu2, tmp, tmp2;
  double *d2min; 
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

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
  which = (int *) R_alloc((size_t) nk, sizeof(int));

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < npoints) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints) maxchunk = npoints;

    for(; i < maxchunk; i++) {

#ifdef SPATSTAT_DEBUG
      Rprintf("\ni=%d\n", i); 
#endif

    /* initialise nn distances and indices */

      d2minK = hu2;
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
	which[k] = -1;
      }

      xi = x[i];
      yi = y[i];

      /* search backward */
      for(left = i - 1; left >= 0; --left)
      {

#ifdef SPATSTAT_DEBUG
	Rprintf("L");
#endif
	dy = yi - y[left];
	dy2 = dy * dy;
	if(dy2 > d2minK)
	  break;

	dx = x[left] - xi;
	d2 =  dx * dx + dy2;
	if (d2 < d2minK) {
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  which[nk1] = left;
	  /* bubble sort */
	  unsorted = TRUE;
	  for(k = nk1; unsorted && k > 0; k--) {
	    k1 = k - 1;
	    if(d2min[k] < d2min[k1]) {
	      /* swap entries */
	      tmp  = d2min[k1];
	      itmp = which[k1];
	      d2min[k1] = d2min[k];
	      which[k1] = which[k];
	      d2min[k] = tmp;
	      which[k] = itmp;
	    } else {
	      unsorted = FALSE;
	    }
	  }
	  /* adjust maximum distance */
	  d2minK = d2min[nk1];
	}
      }

      /* search forward */
      for(right = i + 1; right < npoints; ++right)
	{

#ifdef SPATSTAT_DEBUG
	  Rprintf("R");
#endif
	  dy = y[right] - yi;
	  dy2 = dy * dy;
	  if(dy2 > d2minK)
	    break;

	  dx = x[right] - xi;
	  d2 =  dx * dx + dy2;
	  if (d2 < d2minK) {
	    /* overwrite last entry */
	    d2min[nk1] = d2;
	    which[nk1] = right;
	    /* bubble sort */
	    unsorted = TRUE;
	    for(k = nk1; unsorted && k > 0; k--) {
	      k1 = k - 1;
	      if(d2min[k] < d2min[k1]) {
		/* swap entries */
		tmp  = d2min[k1];
		itmp = which[k1];
		d2min[k1] = d2min[k];
		which[k1] = which[k];
		d2min[k] = tmp;
		which[k] = itmp;
	      } else {
		unsorted = FALSE;
	      }
	    }
	    /* adjust maximum distance */
	    d2minK = d2min[nk1];
	  }
	}

      /* search finished for point i */

#ifdef SPATSTAT_DEBUG
      Rprintf("\n");
#endif

      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
	nnd[nk * i + k] = sqrt(d2min[k]);
	nnwhich[nk * i + k] = which[k] + 1;  /* R indexing */
      }

      /* end of i loop */
    }
  }
}


