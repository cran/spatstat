/*

  nnMDdist.c

  Nearest Neighbour Distances in m dimensions

  $Revision: 1.3 $     $Date: 2011/05/17 12:40:09 $

  Argument x is an m * n matrix 
  with columns corresponding to points
  and rows corresponding to coordinates.

  THE FOLLOWING FUNCTIONS ASSUME THAT THE ROWS OF x 
  ARE SORTED IN ASCENDING ORDER OF THE FIRST COLUMN

  nndMD     Nearest neighbour distances 
  nnwMD     Nearest neighbours and their distances
  nnXwMD    Nearest neighbour from one list to another
  nnXxMD    Nearest neighbour from one list to another, with overlaps

  knndMD    k-th nearest neighbour distances
  knnwMD    k-th nearest neighbours and their distances
*/

#define SPATSTAT_DEBUG 42

#include <R.h>
#include <math.h>
/* #include <stdio.h> */

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

double sqrt();

void nndMD(n, m, x, nnd, huge)
     /* inputs */
     int *n, *m;
     double *x, *huge;
     /* output */
     double *nnd;
{ 
  int npoints, mdimen, i, j, left, right, leftpos, rightpos;
  double dmin, d2, d2min, hu, hu2, xi0, dx0, dxj;
  double *xi;

  npoints = *n;
  mdimen  = *m;
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

#ifdef SPATSTAT_DEBUG
  FILE *out; 
#endif

  hu = *huge;
  hu2 = hu * hu;

#ifdef SPATSTAT_DEBUG
  out = fopen("outnndMD.txt", "w");
#endif

  for(i = 0; i < npoints; i++) {

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\ni=%d\n", i); 
#endif

    dmin = hu;
    d2min = hu2;

    for(j = 0; j < mdimen; j++)
      xi[j] = x[i * mdimen + j];
    xi0 = xi[0];

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n (");
    for(j = 0; j < mdimen; j++)
      fprintf(out, "%lf, ", x[i * mdimen + j]);
    fprintf(out, ")\n");
#endif

    
    /* search backward */
    if(i > 0) {
      for(left = i - 1;
	  left >= 0 && (dx0 = (xi0 - x[left * mdimen])) < dmin ;
	  --left)
	{

#ifdef SPATSTAT_DEBUG
	  fprintf(out, "L=%d, dmin=%lf\n", left, dmin);
#endif

	  d2 = dx0 * dx0;
	  if(mdimen > 1) {
	    leftpos = left * mdimen;
	    for(j = 1; j < mdimen && d2 < d2min; j++) {
	      dxj = xi[j] - x[leftpos + j];
	      d2 += dxj * dxj;
	    }
	  }

	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tupdating dmin=%lf\n", dmin);
#endif
	  }
	}
    }

    /* search forward */
    if(i < npoints - 1) {
      for(right = i + 1;
	  right < npoints && (dx0 = (x[right * mdimen] - xi0)) < dmin ;
	  ++right)
	{

#ifdef SPATSTAT_DEBUG
	  fprintf(out, "R=%d, dmin=%lf\n", right, dmin);
#endif
	  d2 = dx0 * dx0;
	  if(mdimen > 1) {
	    rightpos = right * mdimen;
	    for(j = 1; j < mdimen && d2 < d2min; j++) {
	      dxj = xi[j] - x[rightpos + j];
	      d2 += dxj * dxj;
	    }
	  }

	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tupdating dmin=%lf\n", dmin);
#endif
	  }
	}
    }
#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n");
#endif

    nnd[i] = dmin;
  }

#ifdef SPATSTAT_DEBUG
  fclose(out);
#endif

}

/* nnwMD: same as nndMD, 
   but also returns id of nearest neighbour 
*/

void nnwMD(n, m, x, nnd, nnwhich, huge)
     /* inputs */
     int *n, *m;
     double *x, *huge;
     /* output */
     double *nnd;
     int *nnwhich;
{ 
  int npoints, mdimen, i, j, left, right, leftpos, rightpos, which;
  double dmin, d2, d2min, hu, hu2, xi0, dx0, dxj;
  double *xi;

  npoints = *n;
  mdimen  = *m;
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

#ifdef SPATSTAT_DEBUG
  FILE *out; 
#endif

  hu = *huge;
  hu2 = hu * hu;

#ifdef SPATSTAT_DEBUG
  out = fopen("outnnwMD.txt", "w");
#endif

  for(i = 0; i < npoints; i++) {

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\ni=%d\n", i); 
#endif

    dmin = hu;
    d2min = hu2;
    which = -1;

    for(j = 0; j < mdimen; j++)
      xi[j] = x[i * mdimen + j];
    xi0 = xi[0];

    /* search backward */
    if(i > 0) {
      for(left = i - 1;
	  left >= 0 && (dx0 = (xi0 - x[left * mdimen])) < dmin ;
	  --left)
	{

#ifdef SPATSTAT_DEBUG
	  fprintf(out, "L");
#endif

	  d2 = dx0 * dx0;
	  if(mdimen > 1) {
	    leftpos = left * mdimen;
	    for(j = 1; j < mdimen && d2 < d2min; j++) {
	      dxj = xi[j] - x[leftpos + j];
	      d2 += dxj * dxj;
	    }
	  }

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
	  right < npoints && (dx0 = (x[right * mdimen] - xi0)) < dmin ;
	  ++right)
	{

#ifdef SPATSTAT_DEBUG
	  fprintf(out, "R");
#endif
	  d2 = dx0 * dx0;
	  if(mdimen > 1) {
	    rightpos = right * mdimen;
	    for(j = 1; j < mdimen && d2 < d2min; j++) {
	      dxj = xi[j] - x[rightpos + j];
	      d2 += dxj * dxj;
	    }
	  }

	  if (d2 < d2min) {
	    d2min = d2;
	    dmin = sqrt(d2);
	    which = right;
	  }
	}
    }
#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n");
#endif

    nnd[i] = dmin;
    nnwhich[i] = which;
  }

#ifdef SPATSTAT_DEBUG
  fclose(out);
#endif

}

/* 
   nnXwMD:  for TWO point patterns X and Y,
              find the nearest neighbour 
	      (from each point of X to the nearest point of Y)
	      returning both the distance and the identifier

   Requires both patterns to be sorted in order of increasing z coord
*/

void nnXwMD(m, n1, x1, n2, x2, nnd, nnwhich, huge)
     /* inputs */
     int *m, *n1, *n2;
     double *x1, *x2, *huge;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int mdimen, npoints1, npoints2, i, ell, jleft, jright, jwhich, lastjwhich;
  double dmin, d2, d2min, x1i0, dx0, dxell, hu, hu2;
  double *x1i;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;
  mdimen   = *m;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  x1i = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx  = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  lastjwhich = 0;

  for(i = 0; i < npoints1; i++) {
    dmin = hu;
    d2min = hu2;
    jwhich = -1;
    for(ell = 0; ell < mdimen; ell++) 
      x1i[ell] = x1[i * mdimen + ell];
    x1i0 = x1i[0];

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {
      for(jleft = lastjwhich - 1;
	  jleft >= 0 && (dx0 = (x1i0 - x2[jleft * mdimen])) < dmin ;
	  --jleft)
	{
	  d2 = dx0 * dx0;
	  if(mdimen > 1) {
	    for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
	      dxell = x1i[ell] - x2[jleft * mdimen + ell];
	      d2 += dxell * dxell;
	    }
	  }
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
	  jright < npoints2 && (dx0 = (x2[jright * mdimen] - x1i0)) < dmin ;
	  ++jright)
	{
	  d2 = dx0 * dx0;
	  if(mdimen > 1) {
	    for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
	      dxell = x1i[ell] - x2[jright * mdimen + ell];
	      d2 += dxell * dxell;
	    }
	  }
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
   nnXxMD:  similar to nnXwMD
              but allows X and Y to include common points
	      (which are not to be counted as neighbours)

   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

   Requires both patterns to be sorted in order of increasing y coord
*/

void nnXxMD(m, n1, x1, id1, n2, x2, id2, nnd, nnwhich, huge)
     /* inputs */
     int *m, *n1, *n2;
     double *x1, *x2, *huge;
     int *id1, *id2;
     /* outputs */
     double *nnd;
     int *nnwhich;
{ 
  int mdimen, npoints1, npoints2, i, ell, jleft, jright, jwhich, lastjwhich, id1i;
  double dmin, d2, d2min, x1i0, dx0, dxell, hu, hu2;
  double *x1i;

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;
  mdimen   = *m;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  x1i = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx  = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  lastjwhich = 0;

  for(i = 0; i < npoints1; i++) {
    dmin = hu;
    d2min = hu2;
    jwhich = -1;
    id1i   = id1[i];
    for(ell = 0; ell < mdimen; ell++) 
      x1i[ell] = x1[i * mdimen + ell];
    x1i0 = x1i[0];

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {
      for(jleft = lastjwhich - 1;
	  jleft >= 0 && (dx0 = (x1i0 - x2[jleft * mdimen])) < dmin ;
	  --jleft)
	{
	  /* do not compare identical points */
	  if(id2[jleft] != id1i) {
	    d2 = dx0 * dx0;
	    if(mdimen > 1) {
	      for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
		dxell = x1i[ell] - x2[jleft * mdimen + ell];
		d2 += dxell * dxell;
	      }
	    }
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
	  jright < npoints2 && (dx0 = (x2[jright * mdimen] - x1i0)) < dmin ;
	  ++jright)
	{
	  /* do not compare identical points */
	  if(id2[jright] != id1i) {	  
	    d2 = dx0 * dx0;
	    if(mdimen > 1) {
	      for(ell = 1; ell < mdimen && d2 < d2min; ell++) {
		dxell = x1i[ell] - x2[jright * mdimen + ell];
		d2 += dxell * dxell;
	      }
	    }
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
   knndMD

   nearest neighbours 1:kmax

*/

void knndMD(n, m, kmax, x, nnd, huge)
     /* inputs */
     int *n, *m, *kmax;
     double *x, *huge;
     /* output matrix (kmax * npoints) */
     double *nnd;
{ 
  int npoints, mdimen, nk, nk1, i, j, k, k1, left, right, unsorted;
  double d2, dminK, d2minK, xi0, dx0, dxj, hu, hu2, tmp, tmp2;
  double *dmin, *d2min, *xi;

#ifdef SPATSTAT_DEBUG
  FILE *out; 
#endif

  hu = *huge;
  hu2 = hu * hu;

#ifdef SPATSTAT_DEBUG
  out = fopen("outknndMD.txt", "w");
#endif

  npoints = *n;
  mdimen  = *m;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances
     for the current point
  */

  dmin = (double *) R_alloc((size_t) nk, sizeof(double));
  d2min = (double *) R_alloc((size_t) nk, sizeof(double));

  /* 
     scratch space
  */
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  /* loop over points */

  for(i = 0; i < npoints; i++) {

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\ni=%d\n", i); 
#endif

    /* initialise nn distances */

    dminK  = hu;
    d2minK = hu2;
    for(k = 0; k < nk; k++) {
      dmin[k] = hu;
      d2min[k] = hu2;
    }

    for(j = 0; j < mdimen; j++)
      xi[j] = x[i* mdimen + j];
    xi0 = xi[0];

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n (");
    for(j = 0; j < mdimen; j++)
      fprintf(out, "%lf, ", xi[j]);
    fprintf(out, ")\n");
#endif

    /* search backward */
    for(left = i - 1;
        left >= 0 && (dx0 = (xi0 - x[left * mdimen])) < dminK ;
	--left)
      {

	d2 = dx0 * dx0; 
#ifdef SPATSTAT_DEBUG
	fprintf(out, "L=%d\n", left);
	fprintf(out, "\t 0 ");
#endif
	if(mdimen > 1) {
	  for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	    fprintf(out, "%d ", j);
#endif
	    dxj = xi[j] - x[left * mdimen + j];
	    d2 += dxj * dxj;
	  }
	}
#ifdef SPATSTAT_DEBUG
	fprintf(out, "\n\t sqrt(d2)=%lf\n", sqrt(d2));
#endif
	if (d2 < d2minK) {
	  /* overwrite last entry */
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tsqrt(d2)=%lf overwrites dmin[%d] = %lf\n", 
		  sqrt(d2), nk1, dmin[nk1]);
#endif
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
#endif
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
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
#endif
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dx0 = (x[right * mdimen] - xi0)) < dminK ;
	++right)
      {

#ifdef SPATSTAT_DEBUG
	fprintf(out, "R=%d\n", right);
	fprintf(out, "\t 0 ");
#endif
	d2 = dx0 * dx0; 
	if(mdimen > 1) {
	  for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	    fprintf(out, "%d ", j);
#endif
	    dxj = xi[j] - x[right * mdimen + j];
	    d2 += dxj * dxj;
	  }
	}
#ifdef SPATSTAT_DEBUG
	fprintf(out, "\n\t sqrt(d2)=%lf\n", sqrt(d2));
#endif
	if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tsqrt(d2)=%lf overwrites dmin[%d] = %lf\n", 
		  sqrt(d2), nk1, dmin[nk1]);
#endif
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "sqrt(d2)=%lf\n", sqrt(d2));
#endif
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
#endif
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
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
#endif
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n");
#endif

    /* copy nn distances for point i 
       to output matrix in ROW MAJOR order
    */
    for(k = 0; k < nk; k++) {
      nnd[nk * i + k] = dmin[k];
    }
  }

#ifdef SPATSTAT_DEBUG
  fclose(out);
#endif

}

/* 
   knnwMD

   nearest neighbours 1:kmax

   returns distances and indices

*/


void knnwMD(n, m, kmax, x, nnd, nnwhich, huge)
     /* inputs */
     int *n, *m, *kmax;
     double *x, *huge;
     /* output matrix (kmax * npoints) */
     double *nnd;
     int *nnwhich;
{ 
  int npoints, mdimen, nk, nk1, i, j, k, k1, left, right, unsorted, itmp;
  double d2, dminK, d2minK, xi0, dx0, dxj, hu, hu2, tmp, tmp2;
  double *dmin, *d2min, *xi;
  int *which;

#ifdef SPATSTAT_DEBUG
  FILE *out; 
#endif

  hu = *huge;
  hu2 = hu * hu;

#ifdef SPATSTAT_DEBUG
  out = fopen("outknnwMD.txt", "w");
#endif

  npoints = *n;
  mdimen  = *m;
  nk      = *kmax;
  nk1     = nk - 1;

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  dmin = (double *) R_alloc((size_t) nk, sizeof(double));
  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
  which = (int *) R_alloc((size_t) nk, sizeof(int));

  /* 
     scratch space
  */
  xi = (double *) R_alloc((size_t) mdimen, sizeof(double));
  /*  dx = (double *) R_alloc((size_t) mdimen, sizeof(double)); */

  /* loop over points */

  for(i = 0; i < npoints; i++) {

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\ni=%d\n", i); 
#endif

    /* initialise nn distances */

    dminK  = hu;
    d2minK = hu2;
    for(k = 0; k < nk; k++) {
      dmin[k] = hu;
      d2min[k] = hu2;
      which[k] = -1;
    }

    for(j = 0; j < mdimen; j++)
      xi[j] = x[i* mdimen + j];
    xi0 = xi[0];

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n (");
    for(j = 0; j < mdimen; j++)
      fprintf(out, "%lf, ", x[i * mdimen + j]);
    fprintf(out, ")\n");
#endif

    /* search backward */
    for(left = i - 1;
        left >= 0 && (dx0 = (xi0 - x[left * mdimen])) < dminK ;
	--left)
      {

#ifdef SPATSTAT_DEBUG
	  fprintf(out, "L=%d, dminK=%lf\n", left, dminK);
	  fprintf(out, "\t 0 ");
#endif

	d2 = dx0 * dx0; 
	if(mdimen > 1) {
	  for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	    fprintf(out, "%d ", j);
#endif
	    dxj = xi[j] - x[left * mdimen + j];
	    d2 += dxj * dxj;
	  }
	}
	if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tsqrt(d2)=%lf overwrites dmin[%d] = %lf\n", 
		  sqrt(d2), nk1, dmin[nk1]);
#endif
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  which[nk1] = left;
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
	  fprintf(out, "\twhich[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%d, ", which[k]);
	  fprintf(out, "\n");
#endif
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
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
	  fprintf(out, "\twhich[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%d, ", which[k]);
	  fprintf(out, "\n");
#endif
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

    /* search forward */
    for(right = i + 1;
	right < npoints && (dx0 = (x[right * mdimen] - xi0)) < dminK ;
	++right)
      {

#ifdef SPATSTAT_DEBUG
	  fprintf(out, "R=%d, dminK=%lf\n", right, dminK);
	  fprintf(out, "\t 0 ");
#endif
	d2 = dx0 * dx0; 
	if(mdimen > 1) {
	  for(j = 1; j < mdimen && d2 < d2minK; j++) {
#ifdef SPATSTAT_DEBUG
	    fprintf(out, "%d ", j);
#endif
	    dxj = xi[j] - x[right * mdimen + j];
	    d2 += dxj * dxj;
	  }
	}
	if (d2 < d2minK) {
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tsqrt(d2)=%lf overwrites dmin[%d] = %lf\n", 
		  sqrt(d2), nk1, dmin[nk1]);
#endif
	  /* overwrite last entry */
	  d2min[nk1] = d2;
	  dmin[nk1] = sqrt(d2);
	  which[nk1] = right;
	  /* bubble sort */
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
	  fprintf(out, "\twhich[] before bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%d, ", which[k]);
	  fprintf(out, "\n");
#endif
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
#ifdef SPATSTAT_DEBUG
	  fprintf(out, "\tdmin[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%lf, ", dmin[k]);
	  fprintf(out, "\n");
	  fprintf(out, "\twhich[] after bubble sort:");
	  for(k = 0; k < nk; k++)
	    fprintf(out, "%d, ", which[k]);
	  fprintf(out, "\n");
#endif
	  /* adjust maximum distance */
	  dminK  = dmin[nk1];
	  d2minK = d2min[nk1];
	}
      }

#ifdef SPATSTAT_DEBUG
    fprintf(out, "\n");
#endif

    /* copy nn distances for point i 
       to output matrix in ROW MAJOR order
    */
    for(k = 0; k < nk; k++) {
      nnd[nk * i + k] = dmin[k];
      nnwhich[nk * i + k] = which[k];
    }
  }

#ifdef SPATSTAT_DEBUG
  fclose(out);
#endif

}


