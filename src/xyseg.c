/*

  xyseg.c

  Computation with line segments

  xysegint     compute intersections between line segments

  $Revision: 1.5 $     $Date: 2007/03/16 12:49:51 $

 */

#include <math.h>
#include <stdio.h>

#define NIETS -1.0

/*  
   xysegint

   Determines intersections between each pair of line segments
   drawn from two lists of line segments.

   Line segments are given as x0, y0, dx, dy
   where (x0,y0) is the first endpoint and (dx, dy) is the vector
   from the first to the second endpoint.
   Points along a line segment are represented in parametric
   coordinates, 
            (x,y) = (x0, y0) + t * (dx, dy).

   Output from xysegint() consists of five matrices xx, yy, ta, tb, ok.
   The (i,j)-th entries in these matrices give information about the
   intersection between the i-th segment in list 'a' and the
   j-th segment in list 'b'. The information is

       ok[i,j]  = 1 if there is an intersection
                = 0 if not

       xx[i,j]  = x coordinate of intersection

       yy[i,j]  = y coordinate of intersection

       ta[i,j] = parameter of intersection point
                 relative to i-th segment in list 'a'

       tb[i,j] = parameter of intersection point
                 relative to j-th segment in list 'b'

*/
	     

void xysegint(na, x0a, y0a, dxa, dya, 
              nb, x0b, y0b, dxb, dyb, 
	      eps,
              xx, yy, ta, tb, ok)
     /* inputs (vectors of coordinates) */
     int *na, *nb;
     double *x0a, *y0a, *dxa, *dya, *x0b, *y0b, *dxb, *dyb;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *xx, *yy, *ta, *tb;
     int *ok;
{ 
  int i, j, ma, mb, ijpos;
  double determinant, absdet, diffx, diffy, tta, ttb;

  ma = *na;
  mb = *nb;
  for (j=0; j < mb; j++) 
    {
      for(i = 0; i < ma; i++) {
	ijpos = j * ma + i;
	ok[ijpos] = 0;
	xx[ijpos] = yy[ijpos] = ta[ijpos] = tb[ijpos] = NIETS;
	determinant = dxb[j] * dya[i] - dyb[j] * dxa[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > *eps) {
	  diffx = (x0b[j] - x0a[i])/determinant;
	  diffy = (y0b[j] - y0a[i])/determinant;
	  ta[ijpos] = tta = - dyb[j] * diffx + dxb[j] * diffy;
	  tb[ijpos] = ttb = - dya[i] * diffx + dxa[i] * diffy;
	  if(tta >= 0.0 && tta <= 1.0 && ttb >= 0.0 && ttb <= 1.0) {
	    /* intersection */
	    ok[ijpos] = 1;
	    xx[ijpos] = x0a[i] + tta * dxa[i];
	    yy[ijpos] = y0a[i] + tta * dya[i];
	  }
	}
      }
    }
}
	 
void xysi(na, x0a, y0a, dxa, dya, 
              nb, x0b, y0b, dxb, dyb, 
	      eps,
              ok)
     /* inputs (vectors of coordinates) */
     int *na, *nb;
     double *x0a, *y0a, *dxa, *dya, *x0b, *y0b, *dxb, *dyb;
     /* input (tolerance for determinant) */
     double *eps;  
     /* output (logical) */
     int *ok;
{ 
  int i, j, ma, mb;
  double determinant, absdet, diffx, diffy, tta, ttb;

  *ok = 0;
  ma = *na;
  mb = *nb;
  for (j=0; j < mb; j++) 
    {
      for(i = 0; i < ma; i++) {
	determinant = dxb[j] * dya[i] - dyb[j] * dxa[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > *eps) {
	  diffx = (x0b[j] - x0a[i])/determinant;
	  diffy = (y0b[j] - y0a[i])/determinant;
	  tta = - dyb[j] * diffx + dxb[j] * diffy;
	  ttb = - dya[i] * diffx + dxa[i] * diffy;
	  if(tta >= 0.0 && tta <= 1.0 && ttb >= 0.0 && ttb <= 1.0) {
	    /* intersection */
	    *ok = 1;
	    return;
	  }
	}
      }
    }
}
	 

/*

    Similar to above,
    but computes intersections between all pairs of segments
    in a single list, excluding the diagonal comparisons of course

*/

void xysegXint(n, x0, y0, dx, dy, 
	      eps,
              xx, yy, ti, tj, ok)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *xx, *yy, *ti, *tj;
     int *ok;
{ 
  int i, j, m, mm1, ijpos, jipos, iipos;
  double determinant, absdet, diffx, diffy, tti, ttj;

  m = *n;
 
  mm1 = m - 1;
  for (j=0; j < mm1; j++) 
    {
      for(i = j+1; i < m; i++) {
	ijpos = j * m + i;
	jipos = i * m + j;
	ok[ijpos] = ok[jipos] = 0;
	xx[ijpos] = yy[ijpos] = ti[ijpos] = ti[jipos] = NIETS;
	xx[jipos] = yy[jipos] = tj[ijpos] = tj[jipos] = NIETS;
	determinant = dx[j] * dy[i] - dy[j] * dx[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > *eps) {
	  diffx = (x0[j] - x0[i])/determinant;
	  diffy = (y0[j] - y0[i])/determinant;
	  ti[ijpos] = tti = - dy[j] * diffx + dx[j] * diffy;
	  tj[ijpos] = ttj = - dy[i] * diffx + dx[i] * diffy;
	  tj[jipos] = ti[ijpos];
	  ti[jipos] = tj[ijpos];
	  if(tti >= 0.0 && tti <= 1.0 && ttj >= 0.0 && ttj <= 1.0) {
	    ok[ijpos] = ok[jipos] = 1;
	    xx[ijpos] = xx[jipos] = x0[i] + tti * dx[i];
	    yy[ijpos] = yy[jipos] = y0[i] + tti * dy[i];
	  }
	}
      }
    }

  /* assign diagonal */
  for(i = 0; i < m; i++) {
    iipos = i * m + i;
    ok[iipos] = 0;
    xx[iipos] = yy[iipos] = ti[iipos] = tj[iipos] = NIETS;
  }

}
	 
/*

    Identify self-intersections in a closed polygon

    (Similar to xysegXint,
    but does not compare segments which are cyclically adjacent in the list)

*/

void xypolyselfint(n, x0, y0, dx, dy, 
	      eps,
              xx, yy, ti, tj, ok)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* input (tolerance for determinant) */
     double *eps;  
     /* outputs (matrices) */
     double *xx, *yy, *ti, *tj;
     int *ok;
{ 
  int i, j, m, mm1, mm2, mstop, ijpos, jipos, iipos;
  double determinant, absdet, diffx, diffy, tti, ttj;

  m = *n;

  /* initialise */
  for(i = 0; i < m; i++) 
    for(j = 0; j < m; j++) {
	ijpos = j * m + i;
	ok[ijpos] = 0;
	xx[ijpos] = yy[ijpos] = ti[ijpos] = tj[ijpos] = NIETS;
    }

  if(m <= 2) 
    return;

  /* Compare j with j+2, j+3, ...., m-1
     Don't compare 0 with m-1  */
  mm1 = m - 1;
  mm2 = m - 2;
  for (j=0; j < mm2; j++) 
    {
      mstop = (j > 0) ? m : mm1;
      for(i = j+2; i < mstop; i++) {
	ijpos = j * m + i;
	jipos = i * m + j;
	determinant = dx[j] * dy[i] - dy[j] * dx[i];
	absdet = (determinant > 0) ? determinant : -determinant;
	if(absdet > *eps) {
	  diffx = (x0[j] - x0[i])/determinant;
	  diffy = (y0[j] - y0[i])/determinant;
	  ti[ijpos] = tti = - dy[j] * diffx + dx[j] * diffy;
	  tj[ijpos] = ttj = - dy[i] * diffx + dx[i] * diffy;
	  tj[jipos] = ti[ijpos];
	  ti[jipos] = tj[ijpos];
	  if(tti >= 0.0 && tti <= 1.0 && ttj >= 0.0 && ttj <= 1.0) {
	    ok[ijpos] = ok[jipos] = 1;
	    xx[ijpos] = xx[jipos] = x0[i] + tti * dx[i];
	    yy[ijpos] = yy[jipos] = y0[i] + tti * dy[i];
	  }
	}
      }
    }
}
	 

/*
  Just determines whether there is self-intersection
  (exits quicker & uses less space)
*/


void xypsi(n, x0, y0, dx, dy, xsep, ysep, eps, proper, answer)
     /* inputs (vectors of coordinates) */
     int *n;
     double *x0, *y0, *dx, *dy;
     /* inputs (distances beyond which intersection is impossible) */
     double *xsep, *ysep;
     /* input (tolerance for determinant) */
     double *eps;  
     /* input (flag) */
     int *proper;
     /* output */
     int *answer;
{ 
  int i, j, m, mm1, mm2, mstop, prop;
  double determinant, absdet, diffx, diffy, tti, ttj;
  double Xsep, Ysep;

  m = *n;
  prop = *proper;
  Xsep = *xsep;
  Ysep = *ysep;

  *answer = 0;

  if(m <= 2) 
    return;

  /* Compare j with j+2, j+3, ...., m-1
     Don't compare 0 with m-1  */
  mm1 = m - 1;
  mm2 = m - 2;
  for (j=0; j < mm2; j++) 
    {
      mstop = (j > 0) ? m : mm1;
      for(i = j+2; i < mstop; i++) {
	diffx = x0[j] - x0[i];
	diffy = y0[j] - y0[i];
	if(diffx < Xsep && diffx > -Xsep && diffy < Ysep & diffy > -Ysep) {
	  determinant = dx[j] * dy[i] - dy[j] * dx[i];
	  absdet = (determinant > 0) ? determinant : -determinant;
	  if(absdet > *eps) {
	    diffx = diffx/determinant;
	    diffy = diffy/determinant;
	    tti = - dy[j] * diffx + dx[j] * diffy;
	    ttj = - dy[i] * diffx + dx[i] * diffy;
	    if(tti >= 0.0 && tti <= 1.0 && ttj >= 0.0 && ttj <= 1.0) {
              /* intersection occurs */
	      if(prop == 0 ||
		 (tti != 0.0 && tti != 1.0) || 
		 (ttj != 0.0 && ttj != 1.0)) {
              /* proper intersection */
		*answer = 1;
		return;
	      }
	    }
	  }
	}
      }
    }
}

	 
