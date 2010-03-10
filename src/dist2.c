# include <math.h>

#define YES (0 == 0)
#define NO (!YES)

/* 

   dist2:   squared distance in torus

   dist2thresh: faster code for testing whether dist2 < r2

*/

double dist2(u,v,x,y,period)
     double u, v, x, y, *period;
{
  double wide, high, a1, a2, a3, b1, b2, b3, amin, bmin, d2;

  wide = period[0];
  high = period[1];

  /*

  The following line is commented out, 
  as we now avoid calling 'dist2' in this case:

  if(wide < 0.) return pow(u-x,2) + pow(v-y,2);

  */

  a1 = pow(u-x,2);
  a2 = pow(u-x+wide,2);
  a3 = pow(u-x-wide,2);

  b1 = pow(v-y,2);
  b2 = pow(v-y+high,2);
  b3 = pow(v-y-high,2);

  amin = a1;
  if(a2 < amin) amin=a2;
  if(a3 < amin) amin=a3;

  bmin = b1;
  if(b2 < bmin) bmin=b2;
  if(b3 < bmin) bmin=b3;

  d2 = amin + bmin;
  
  return d2;
}

double dist2either(u,v,x,y,period)
     double u, v, x, y, *period;
{
  if(period[0] < 0.) return pow(u-x,2) + pow(v-y,2);
  return(dist2(u,v,x,y,period));
}

int dist2thresh(u,v,x,y,period,r2)
     double u, v, x, y, *period, r2;
{
  double wide, high, amin, bmin, aa, bb, residue, dx, dy;

  wide = period[0];
  high = period[1];
  
  dx = u - x;
  amin = dx * dx;
  aa = (dx - wide) * (dx - wide);
  if(aa < amin) amin = aa;
  aa = (dx + wide) * (dx + wide);
  if(aa < amin) amin = aa;

  residue = r2 - amin;
  if(residue < 0)
    return NO;

  dy = v - y;
  bmin = dy * dy;
  bb = (dy - high) * (dy - high);
  if(bb < bmin) bmin= bb;
  bb = (dy + high) * (dy + high);
  if(bb < bmin) bmin= bb;

  residue = residue - bmin;

  return (residue >= 0);
}
