/*
       connected.c

       Connected component transform of a discrete binary image
       (8-connected topology)
       
       $Revision: 1.5 $ $Date: 2012/03/27 01:42:56 $

       
*/

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>

#include "raster.h"

void   shape_raster();

void
comcommer(im)
     Raster  *im;            
     /* raster must have been dimensioned by shape_raster() */
     /* Pixel values assumed to be 0 in background, and 
        distinct nonzero integers in foreground */
{
  int	j,k;
  int rmin, rmax, cmin, cmax;
  int label, curlabel, minlabel;
  int nchanged;

  /* image boundaries */
  rmin = im->rmin;
  rmax = im->rmax;
  cmin = im->cmin;
  cmax = im->cmax;

#define ENTRY(ROW, COL) Entry(*im, ROW, COL, int)

#define UPDATE(ROW,COL,BEST,NEW) \
     NEW = ENTRY(ROW, COL); \
     if(NEW != 0 && NEW < BEST) \
       BEST = NEW

  nchanged = 1;

  while(nchanged >0) {
    nchanged = 0;
    R_CheckUserInterrupt();
    for(j = rmin; j <= rmax; j++) {
      for(k = cmin; k <= cmax; k++) {
	curlabel = ENTRY(j, k);
	if(curlabel != 0) {
	  minlabel = curlabel;
	  UPDATE(j-1, k-1, minlabel, label);
	  UPDATE(j-1, k,   minlabel, label);
	  UPDATE(j-1, k+1, minlabel, label);
	  UPDATE(j,   k-1, minlabel, label);
	  UPDATE(j,   k,   minlabel, label);
	  UPDATE(j,   k+1, minlabel, label);
	  UPDATE(j+1, k-1, minlabel, label);
	  UPDATE(j+1, k,   minlabel, label);
	  UPDATE(j+1, k+1, minlabel, label);
	  if(minlabel < curlabel) {
	    ENTRY(j, k) = minlabel;
	    nchanged++;
	  }
	}
      }
    }
  }
}

/* R interface */

void concom(mat, nr, nc)
     int   *mat;        /* input:  binary image */
     int *nr, *nc;      /* raster dimensions
			   EXCLUDING margin of 1 on each side */
{
  Raster im;

  shape_raster( &im, (void *) mat, 
		(double) 1, (double) 1,
		(double) *nc, (double) *nr, 
		*nr+2, *nc+2, 1, 1);
  comcommer(&im);
}	

