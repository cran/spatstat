/*
       exactPdist.c

       `Pseudoexact' distance transform of a discrete binary image
       (the closest counterpart to `exactdist.c')
       
       $Revision: 1.5 $ $Date: 2000/07/11 10:55:13 $

       
*/

#include <math.h>
#include "raster.h"

void   dist_to_bdry();

void
ps_exact_dt(in, dist, row, col)
        Raster  *in;            /* input:  binary image */
	Raster	*dist;		/* output: exact distance to nearest point */
	Raster	*row;		/* output: row index of closest point */
	Raster	*col;		/* output: column index of closest point */
	/* rasters must have been dimensioned by shape_raster()
	   and must all have identical dimensions and margins */
{
	long	i,j,k,l,m;
	double	d, x, y;
	long	r, c;
	double	dnew;
	double  bdiag;
	double  huge;
	long  *ip; 
	double *dp;
	
	    /* initialise */
#define UNDEFINED -1
#define Is_Defined(I) (I >= 0)
#define Is_Undefined(I) (I < 0)
	
	Clear(*row,long,UNDEFINED)
	Clear(*col,long,UNDEFINED)
		
	huge = 2.0 * DistanceSquared(dist->xmin,dist->ymin,dist->xmax,dist->ymax); 
	Clear(*dist,double,huge)


	  /* if input pixel is TRUE, set distance to 0 and make pixel point to itself */
	for(j = in->rmin; j <= in->rmax; j++)
	for(k = in->cmin; k <= in->cmax; k++) 
	  if(Entry(*in, j, k, long) != 0) {
	      Entry(*dist, j, k, double) = 0.0;
	      Entry(*row,  j, k, long)   = j;
	      Entry(*col,  j, k, long)   = k;
	  }

	/* how to update the distance values */
	
#define GETVALUES(ROW,COL) \
	x = Xpos(*in, COL); \
	y = Ypos(*in, ROW); \
	d = Entry(*dist,ROW,COL,double); 

#define COMPARE(ROW,COL,RR,CC,BOUND) \
	r = Entry(*row,RR,CC,long); \
	c = Entry(*col,RR,CC,long); \
	if(Is_Defined(r) && Is_Defined(c) \
	   && Entry(*dist,RR,CC,double) < d) { \
	     dnew = DistanceSquared(x, y, Xpos(*in,c), Ypos(*in,r)); \
	     if(dnew < d) { \
		Entry(*row,ROW,COL,long) = r; \
		Entry(*col,ROW,COL,long) = c; \
		Entry(*dist,ROW,COL,double) = dnew; \
		d = dnew; \
	     } \
	}

	/* bound on diagonal step distance squared */
	bdiag = (in->xstep * in->xstep + in->ystep * in->ystep);
	
	/* forward pass */

	for(j = in->rmin; j <= in->rmax; j++)
	for(k = in->cmin; k <= in->cmax; k++) {
	        GETVALUES(j, k)
		COMPARE(j,k, j-1,k-1, bdiag)
		COMPARE(j,k, j-1,  k, in->ystep)
		COMPARE(j,k, j-1,k+1, bdiag)
		COMPARE(j,k, j,  k-1, in->xstep)
		  }

	/* backward pass */

	for(j = in->rmax; j >= in->rmin; j--) 
	for(k = in->cmax; k >= in->cmin; k--) {
	        GETVALUES(j, k)
		COMPARE(j,k, j+1,k+1, bdiag)
		COMPARE(j,k, j+1,  k, in->ystep)
		COMPARE(j,k, j+1,k-1, bdiag)
		COMPARE(j,k, j,  k+1, in->xstep)
		  } 

	/* take square roots of distances^2 */

	for(j = in->rmax; j >= in->rmin; j--) 
	for(k = in->cmax; k >= in->cmin; k--) 
	        Entry(*dist,j,k,double) = sqrt(Entry(*dist,j,k,double));

}

/* S interface */

ps_exact_dt_S(xmin, ymin, xmax, ymax, nr, nc,
	   in, distances, rows, cols, boundary)
	double *xmin, *ymin, *xmax, *ymax;  	  /* x, y dimensions */
	long *nr, *nc;	 	                  /* raster dimensions
				                     EXCLUDING margin of 1 on each side */
	long   *in;              /* input:  binary image */
	double *distances;	/* output: distance to nearest point */
	long   *rows;	        /* output: row of nearest point (start= 0) */
	long   *cols;	        /* output: column of nearest point (start = 0) */
	double *boundary;       /* output: distance to boundary of rectangle */
	/* all images must have identical dimensions including a margin of 1 on each side */
{
	Raster data, dist, row, col, bdist;

	shape_raster( &data, (char *) in, *xmin,*ymin,*xmax,*ymax,
			    *nr+2, *nc+2, 1, 1);
	shape_raster( &dist, (char *) distances,*xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	shape_raster( &row, (char *) rows, *xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	shape_raster( &col, (char *) cols, *xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	shape_raster( &bdist, (char *) boundary, *xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	
	ps_exact_dt(&data, &dist, &row, &col);

	dist_to_bdry(&bdist);
}	
