/*
      raster.h

      Definition of raster structures & operations

      requires <math.h> (for floor())

      $Revision: 1.2 $ $Date: 1997/08/18 16:03:31 $
*/

typedef struct Raster{
 /* array of data */
	char		*data;		/* coerced to appropriate type */
	long	nrow;		/* dimensions of entire array */
	long	ncol;
	unsigned long	length;
	long	rmin;		/* position of valid subrectangle */
	long	rmax;
	long	cmin;
	long	cmax;
/* definition of mapping into continuous space */
	double	x0;	/* position of entry (rmin,cmin) */
	double	y0;
	double	x1;	/* position of entry (rmax,cmax) */
	double	y1;
	double	xstep;	/* x increment for each column step */
	double	ystep;	/* y increment for each row step */
	                /*
			   xstep = (x1 - x0)/(cmax - cmin)
			         = (x1 - x0)/(number of valid columns - 1)
			   CAN BE POSITIVE OR NEGATIVE 
			 */
	 /* image of valid subrectangle */
	double	xmin;	/* = min{x0,x1} */
	double	xmax;
	double	ymin;
	double	ymax;
} Raster;

/*      how to clear the data      */

#define Clear(ARRAY,TYPE,VALUE) \
       { unsigned long i; TYPE *p; \
	 for(i = 0, p = (TYPE *) (ARRAY).data; i < (ARRAY).length; i++, p++) \
	 *p = VALUE; }
		
/* 	how to index a rectangular array
	stored sequentially in row-major order */

#define Entry(ARRAY,ROW,COL,TYPE) \
	((TYPE *)((ARRAY).data))[COL + (ROW) * ((ARRAY).ncol)]

     /* test for indices inside subrectangle */
	
#define Inside(ARRAY,ROW,COL) \
	( (ROW >= (ARRAY).rmin) && (ROW <= (ARRAY).rmax) && \
	(COL >= (ARRAY).cmin) && (COL <= (ARRAY).cmax))

     /* how to compute the position in R^2 corresponding to a raster entry */

#define Xpos(ARRAY,COL) \
	((ARRAY).x0 + (ARRAY).xstep * (COL - (ARRAY).cmin))
#define Ypos(ARRAY,ROW) \
	((ARRAY).y0 + (ARRAY).ystep * (ROW - (ARRAY).rmin))

#define Distance(X,Y,XX,YY) sqrt((X - XX)* (X - XX) + (Y - YY) * (Y - YY))

#define DistanceTo(X,Y,ARRAY,ROW,COL)\
	Distance(X,Y,Xpos(ARRAY,COL),Ypos(ARRAY,ROW))

#define DistanceSquared(X,Y,XX,YY) ((X - XX)* (X - XX) + (Y - YY) * (Y - YY))

#define DistanceToSquared(X,Y,ARRAY,ROW,COL)\
	DistanceSquared(X,Y,Xpos(ARRAY,COL),Ypos(ARRAY,ROW))


  /* how to map a point (x,y) in R^2 to a raster entry */
  /*
     (x,y) is guaranteed to lie in the rectangle bounded by
     the images of the entries (r,c), (r+1,c), (r,c+1), (r+1,c+1)
     where r = RowIndex(..) and c = ColIndex(..).
  */

#define RowIndex(ARRAY,Y) \
	((ARRAY).rmin + (long) floor(((Y) - (ARRAY).y0)/(ARRAY).ystep))
#define ColIndex(ARRAY,X) \
	((ARRAY).cmin + (long) floor(((X) - (ARRAY).x0)/(ARRAY).xstep))
	
