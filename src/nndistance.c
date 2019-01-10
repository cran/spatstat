/*

  nndistance.c

  Nearest Neighbour Distances between points

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2012
  Licence: GNU Public Licence >= 2

  $Revision: 1.22 $     $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  SINGLE LIST:
  nndistsort    Nearest neighbour distances 
  nnwhichsort   Nearest neighbours
  nnsort        Nearest neighbours & distances

  ONE LIST TO ANOTHER LIST:
  nnXdist       Nearest neighbour distance from one list to another
  nnXwhich      Nearest neighbour ID from one list to another
  nnX           Nearest neighbour ID & distance from one list to another

  ONE LIST TO ANOTHER OVERLAPPING LIST:
  nnXEdist      Nearest neighbour distance from one list to another, overlapping
  nnXEwhich     Nearest neighbour ID from one list to another, overlapping
  nnXE          Nearest neighbour ID & distance 

*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

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

/* general interface */

void nnXinterface(n1, x1, y1, id1, 
		  n2, x2, y2, id2, 
		  exclude, wantdist, wantwhich,
		  nnd, nnwhich, 
		  huge)
     /* inputs */
     int *n1, *n2;
     double *x1, *y1, *x2, *y2, *huge;
     int *id1, *id2;
     /* options */
     int *exclude, *wantdist, *wantwhich;
     /* outputs */
     double *nnd;
     int *nnwhich;
{
  void nnX(), nnXdist(), nnXwhich();
  void nnXE(), nnXEdist(), nnXEwhich();
  int ex, di, wh;
  ex = (*exclude != 0);
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  if(!ex) {
    if(di && wh) {
      nnX(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } else if(di) {
      nnXdist(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } else if(wh) {
      nnXwhich(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } 
  } else {
    if(di && wh) {
      nnXE(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } else if(di) {
      nnXEdist(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } else if(wh) {
      nnXEwhich(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } 
  }
}


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

