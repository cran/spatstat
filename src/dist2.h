/* 
   dist2.h 

   External declarations for the functions defined in dist2.c
   and
   In-line cpp macros for similar purposes

   $Revision: 1.10 $ $Date: 2013/01/14 11:25:49 $

*/

double dist2(double u, double v, double x, double y, double *period);

double dist2either(double u, double v, double x, double y, double *period);

int dist2thresh(double u, double v, double x, double y, double *period, double r2);

int dist2Mthresh(double u, double v, double x, double y, double *period, double r2);

/* 
    IF_CLOSE tests whether (U,V) and (X,Y) are closer than squared distance R2.

    IF_CLOSE_D2 does the same test, but also computes the squared distance D2

*/

#define IF_CLOSE(U,V,X,Y,R2) \
  X_IF_CLOSE(U,V,X,Y,R2, scratchDX, scratchDY, scratchRESIDUE)	

#define X_IF_CLOSE(U,V,X,Y,R2,DX,DY,RESIDUE)	\
  { register double DX, DY, RESIDUE;  \
    DX = X - U;	\
    RESIDUE = R2 - DX * DX;	   \
    if(RESIDUE >= 0.0) {	   \
    DY = Y - V; \
    if(RESIDUE >= DY * DY) 

#define END_IF_CLOSE }}

#define IF_CLOSE_D2(U,V,X,Y,R2,D2) \
  X_IF_CLOSE_D2(U,V,X,Y,R2,D2, scratchDX, scratchDY, scratchDX2)	

#define X_IF_CLOSE_D2(U,V,X,Y,R2,D2,DX,DY,DX2)	\
  { register double DX, DY, DX2;  \
    DX = X - U;	\
    DX2 = DX * DX;	   \
    if(DX2 <= R2) {	   \
    DY = Y - V; \
    D2 = DX2 + DY * DY; \
    if(D2 <= R2) 

#define END_IF_CLOSE_D2 }}

/* 
   counterparts for periodic distance
*/

#define IF_CLOSE_PERIODIC(U,V,X,Y,PERIOD,R2) \
  X_IF_CLOSE_PERIODIC(U,V,X,Y,PERIOD,R2,\
		      scratchDX,scratchDY,scratchDXP,scratchDYP,scratchRESIDUE)

#define X_IF_CLOSE_PERIODIC(U,V,X,Y,PERIOD,R2,DX,DY,DXP,DYP,RESIDUE)	\
  { register double DX, DY, DXP, DYP, RESIDUE;			\
    DX = X - U;	\
    if(DX < 0.0) DX = -DX; \
    DXP = (PERIOD)[0] - DX; \
    DX = (DX < DXP) ? DX : DXP; \
    RESIDUE = R2 - DX * DX;    \
    if(RESIDUE >= 0.0) {       \
      DY = V - Y;              \
      if(DY < 0.0) DY = -DY;    \
      DYP = (PERIOD)[1] - DY;   \
      DY = (DY < DYP) ? DY : DYP;  \
      if (RESIDUE >= DY * DY) 

#define END_IF_CLOSE_PERIODIC }}

#define IF_CLOSE_PERIODIC_D2(U,V,X,Y,PERIOD,R2,D2)			\
  X_IF_CLOSE_PERIODIC_D2(U,V,X,Y,PERIOD,R2,D2,\
			 scratchDX,scratchDY, scratchDXP, scratchDYP)

#define X_IF_CLOSE_PERIODIC_D2(U,V,X,Y,PERIOD,R2,D2,DX,DY,DXP,DYP)	\
  { register double DX, DY, DXP, DYP;				\
    DX = X - U;	\
    if(DX < 0.0) DX = -DX; \
    DXP = (PERIOD)[0] - DX; \
    DX = (DX < DXP) ? DX : DXP; \
    D2 = DX * DX;    \
    if(D2 <= R2) {       \
      DY = V - Y;              \
      if(DY < 0.0) DY = -DY;    \
      DYP = (PERIOD)[1] - DY;   \
      DY = (DY < DYP) ? DY : DYP;  \
      D2 += DY * DY; \
      if (D2 <= R2) 

#define END_IF_CLOSE_PERIODIC_D2 }}

