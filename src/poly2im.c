/*
  poly2im.c

  Conversion from (x,y) polygon to pixel image

  $Revision$ $Date$

*/

#include <math.h>

#define OUT(I,J) out[I + (J) * Ny]

void 
poly2im(xp, yp, np, nx, ny, out) 
     double *xp, *yp; /* polygon vertices, anticlockwise, CLOSED  */
     int *np; 
     int *nx, *ny; /* INTEGER raster points from (0,0) to (nx-1, ny-1) */
     int *out;  /* output matrix [ny, nx], initialised to zero */
{
  int Np, Nx, Ny, Nxy;
  int i, j, k;
  double x0, y0, x1, y1, xleft, xright, yleft, yright;
  double dx, dy, y, slope, intercept;
  int jleft, jright, imax;
  int sign;

  Np = *np;
  Nx = *nx;
  Ny = *ny;
  Nxy = Nx * Ny;

  /* run through polygon edges */
  for(k = 0; k < Np-1; k++) {
    x0 = xp[k];
    y0 = yp[k];
    x1 = xp[k+1];
    y1 = yp[k+1];
    if(x0 < x1) {
      xleft = x0;
      xright = x1;
      yleft = y0;
      yright = y1;
      sign = -1;
    } else {
      xleft = x1;
      xright = x0;
      yleft = y1;
      yright = y0;
      sign = +1;
    }
    /* determine relevant columns of pixels */
    jleft = (int) ceil(xleft);
    jright = (int) floor(xright);
    if(jleft < Nx && jright >= 0 && jleft <= jright) {
      if(jleft < 0) { jleft = 0; } 
      if(jright >= Nx) {jright = Nx - 1; }
      /* equation of edge */
      dx = xright - xleft;
      dy = yright - yleft;
      slope = dy/dx;
      intercept = yleft - slope * xleft;
      /* visit relevant columns */
      for(j = jleft; j <= jright; j++) {
	y = slope * ((double) j) + intercept;
	imax = (int) floor(y);
	if(imax >= Ny) imax = Ny-1;
	if(imax >= 0) {
	  /* increment entries below edge in this column */
	  for(i = 0; i <= imax; i++) {
	    OUT(i,j) += sign;
	  }
	}
      }
    }
  }
}
