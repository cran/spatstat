/* 

   utils.c

   $Revision: 1.1 $  $Date: 2006/09/14 03:46:27 $

   Small utilities

*/

void drevcumsum(double *x, int *nx) {
  int i;
  double sumx;
  double *xp;
  
  i = *nx - 1;
  xp = x + i;
  sumx = *xp;
  while(i > 0) {
    --i;
    --xp;
    sumx += *xp;
    *xp = sumx;
  }
}

void irevcumsum(int *x, int *nx) {
  int i;
  int sumx;
  int *xp;
  
  i = *nx - 1;
  xp = x + i;
  sumx = *xp;
  while(i > 0) {
    --i;
    --xp;
    sumx += *xp;
    *xp = sumx;
  }
}
