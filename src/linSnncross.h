/* 
   linSnncross.h

   Function body definitions with macros

   Sparse representation of network

   $Revision: 1.5 $  $Date: 2018/12/18 02:43:11 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

   Macros used:
   FNAME   name of function
   WHICH   whether 'nnwhich' is required
   HUH     debugging

   ! Data points must be ordered by segment index !
*/

void 
FNAME(np, sp, tp,  /* data points 'from' (ordered by sp) */
      nq, sq, tq, /* data points 'to'   (ordered by sq) */
      nv, /* number of network vertices */
      ns, from, to,  /* segments */
      seglen,  /* segment lengths */
      huge, /* value taken as infinity */
      tol, /* tolerance for updating distances */
      /* OUTPUT */
#ifdef WHICH
      nndist,  /* nearest neighbour distance for each point */
      nnwhich  /* identifies nearest neighbour */
#else 
      nndist  /* nearest neighbour distance for each point */
#endif
)
  int *np, *nq, *nv, *ns;
  int *from, *to, *sp, *sq; /* integer vectors (mappings) */
  double *tp, *tq; /* fractional location coordinates */
  double *huge, *tol;
  double *seglen; 
  double *nndist; /* nearest neighbour distance for each point */
#ifdef WHICH
  int *nnwhich; /* identifies nearest neighbour */
#endif
{
  int Np, Nq, Nv, i, j, ivleft, ivright, jfirst, jlast, k;
  double d, hugevalue, slen, tpi;
  double *dminvert;  /* min dist from each vertex */
#ifdef WHICH
  int *whichvert;   /* which min from each vertex */
#endif 

  Np = *np;
  Nq = *nq;
  Nv = *nv;
  hugevalue = *huge;

  /* First compute min distance to target set from each vertex */
  dminvert = (double *) R_alloc(Nv, sizeof(double));
#ifdef WHICH
  whichvert = (int *) R_alloc(Nv, sizeof(int));
  Clinvwhichdist(nq, sq, tq, nv, ns, from, to, seglen, huge, tol, 
		 dminvert, whichvert);
#else
  Clinvdist(nq, sq, tq, nv, ns, from, to, seglen, huge, tol, 
	    dminvert);
#endif

#ifdef HUH
  Rprintf("Initialise answer\n");
#endif
  /* initialise nn distances from source points */
  for(i = 0; i < Np; i++) {
    nndist[i] = hugevalue;
#ifdef WHICH
    nnwhich[i] = -1;
#endif
  }

  /* run through all source points */
#ifdef HUH
  Rprintf("Run through source points\n");
#endif
  jfirst = 0;
  for(i = 0; i < Np; i++) {
    tpi = tp[i];
    k = sp[i];   /* segment containing this point */
    slen = seglen[k];
    ivleft = from[k];
    ivright = to[k];
#ifdef HUH
    Rprintf("Source point %d lies on segment %d = [%d,%d]\n", 
	    i, k, ivleft, ivright);
#endif
    d = slen * tpi + dminvert[ivleft];
    if(nndist[i] > d) {
#ifdef HUH
      Rprintf("\tMapping to left endpoint %d, distance %lf\n", ivleft, d);
#endif
      nndist[i] = d;
#ifdef WHICH
      nnwhich[i] = whichvert[ivleft];
#endif
    }
    d = slen * (1.0 - tpi) + dminvert[ivright];
    if(nndist[i] > d) {
#ifdef HUH
      Rprintf("\tMapping to right endpoint %d, distance %lf\n", ivright, d);
#endif
      nndist[i] = d;
#ifdef WHICH
      nnwhich[i] = whichvert[ivright];
#endif
    }
    /* find any target points in this segment */
    while(jfirst < Nq && sq[jfirst] < k) jfirst++;
    jlast = jfirst;
    while(jlast < Nq && sq[jlast] == k) jlast++;
    --jlast;
    /* if there are no such points, then jlast < jfirst */
    if(jfirst <= jlast) {
      for(j = jfirst; j <= jlast; j++) {
	d = slen * fabs(tq[j] - tpi);
	if(nndist[i] > d) {
	  nndist[i] = d;
#ifdef WHICH
	  nnwhich[i] = j;
#endif
	}
      }
    }
  }
}

