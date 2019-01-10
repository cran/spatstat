#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

/* 
  Knone.c

  Efficient computation of uncorrected estimates of K 
  for large datasets
  
  KnoneI()   Estimates K function, 
               returns integer numerator

  KnoneD()   Estimates K function, 
               returns double precision numerator

  Kwnone()   Estimates Kinhom, 
               returns double precision numerator

  Functions require (x,y) data to be sorted in ascending order of x
  and expect r values to be equally spaced and starting at zero
   
  $Revision: 1.3 $ $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

#undef WEIGHTED

#define FNAME KnoneI
#define OUTTYPE int
#include "Knone.h"

#undef FNAME
#undef OUTTYPE
#define FNAME KnoneD
#define OUTTYPE double
#include "Knone.h"

#undef FNAME
#undef OUTTYPE
#define FNAME Kwnone
#define WEIGHTED
#define OUTTYPE double
#include "Knone.h"



