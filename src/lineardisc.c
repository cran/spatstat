#include <R.h>

/* 
   lineardisc.c

   Disc of radius r in linear network

   $Revision: 1.5 $  $Date: 2011/07/26 09:21:18 $

*/

#define DPATH(I,J) dpath[(J) + Nv * (I)]

#define TRUE (0 == 0)
#define FALSE (!TRUE)

void 
lineardisc(f, seg, /* centre of disc (local coords) */
	   r,      /* radius of disc */
	   nv, xv, yv,   /* network vertices */
	   ns, from, to,  /* segments */
	   dpath,  /* shortest path distances between vertices */
	   lengths, /* segment lengths */
	   allinside, boundary, dxv, nendpoints)
     int *nv, *ns;
     int *from, *to; /* integer vectors (mappings) */
     double *f, *r; 
     int *seg;
     double *xv, *yv; /* vectors of coordinates of vertices */
     double *dpath; /* matrix of shortest path distances between vertices */
     double *lengths; /* vector of segment lengths */
     /* OUTPUTS */
     int *allinside, *boundary; /* vectors of status for each segment */
     double *dxv; /* vector of distances for each vertex */
     int *nendpoints;
{
  int Nv, Ns;
  double f0, rad;
  int seg0;

  int i, A, B, fromi, toi, allin, bdry, reachable, nends;
  double length0, dxA, dxB, dxAvi, dxBvi, residue;
  double *resid; 
  int *covered;

  Nv = *nv;
  Ns = *ns;

  f0 = *f;
  seg0 = *seg;
  rad = *r;

  /* endpoints of segment containing centre */
  A = from[seg0];
  B = to[seg0];

  /* distances from x to  A and B */
  length0 = lengths[seg0];
  dxA = f0 * length0;
  dxB = (1-f0) * length0;

  /* visit vertices */
  covered = (int *) R_alloc((size_t) Nv, sizeof(int));
  resid = (double *) R_alloc((size_t) Nv, sizeof(double));
  for(i = 0; i < Nv; i++) {
    /* distance going through A */
    dxAvi = dxA + DPATH(A,i);
    /* distance going through B */
    dxBvi = dxB + DPATH(B,i);
    /* shortest path distance to this vertex */
    dxv[i] = (dxAvi < dxBvi) ? dxAvi : dxBvi;
    /* distance left to 'spend' from this vertex */
    residue = rad - dxv[i];
    resid[i] = (residue > 0)? residue : 0;
    /* determine whether vertex i is inside the disc of radius r */
    covered[i] = (residue >= 0);
  }

  /* 
     Now visit line segments. 
  */
  nends = 0;

  for(i = 0; i < Ns; i++) {
    /* 
     Determine which line segments are completely inside the disc,
     and which cross the boundary.
    */
    if(i == seg0) {
      /* initial segment: disc starts from centre (x, y) */
      allin = covered[A] && covered[B];
      bdry  = !allin;
      if(bdry) {
	if(!covered[A]) nends++;
	if(!covered[B]) nends++;
      }
    } else {
      /* another segment: disc extends in from either endpoint */
      fromi = from[i];
      toi   = to[i];
      reachable = (covered[fromi] || covered[toi]);
      if(reachable) {
	allin = (resid[fromi] + resid[toi] >= lengths[i]);
	bdry = !allin;
      } else allin = bdry = FALSE;
      if(bdry) {
	if(covered[fromi]) nends++;
	if(covered[toi]) nends++;
      }
    }
    allinside[i] = allin;
    boundary[i] = bdry;
  }

  *nendpoints = nends;
}

/* ------------------------------------------------- */
/*   count endpoints of several discs in a network   */
/* ------------------------------------------------- */

void 
countends(np, f, seg, /* centres of discs (local coords) */
	  r,                /* radii of discs */
	  nv, xv, yv,   /* network vertices */
	  ns, from, to,  /* network segments */
	  dpath,  /* shortest path distances between vertices */
	  lengths, /* segment lengths */
	  nendpoints /* output counts of endpoints */
	  )
     int *np, *nv, *ns;
     int *from, *to; /* integer vectors (mappings) */
     double *f, *r; 
     int *seg;
     double *xv, *yv; /* vectors of coordinates of vertices */
     double *dpath; /* matrix of shortest path distances between vertices */
     double *lengths; /* vector of segment lengths */
     /* OUTPUT */
     int *nendpoints;
{
  int Np, Nv, Ns;
  double f0, rad;
  int seg0;

  int i, m, A, B, fromi, toi, allin, bdry, reachable, nends;
  double length0, dxA, dxB, dxAvi, dxBvi, dxvi, residue;
  double *resid; 
  int *covered;

  Np = *np;
  Nv = *nv;
  Ns = *ns;

  covered = (int *) R_alloc((size_t) Nv, sizeof(int));
  resid = (double *) R_alloc((size_t) Nv, sizeof(double));

  /* loop over centre points */
  for(m = 0; m < Np; m++) {
    f0 = f[m];
    seg0 = seg[m];
    rad = r[m];

    /* endpoints of segment containing centre */
    A = from[seg0];
    B = to[seg0];

    /* distances from centre to A and B */
    length0 = lengths[seg0];
    dxA = f0 * length0;
    dxB = (1-f0) * length0;

    /* visit vertices */
    for(i = 0; i < Nv; i++) {
      /* distance going through A */
      dxAvi = dxA + DPATH(A,i);
      /* distance going through B */
      dxBvi = dxB + DPATH(B,i);
      /* shortest path distance to this vertex */
      dxvi = (dxAvi < dxBvi) ? dxAvi : dxBvi;
      /* distance left to 'spend' from this vertex */
      residue = rad - dxvi;
      resid[i] = (residue > 0)? residue : 0;
      /* determine whether vertex i is inside the disc of radius r */
      covered[i] = (residue >= 0);
    }

    /* 
       Now visit line segments. 
    */
    nends = 0;

    for(i = 0; i < Ns; i++) {
      /* 
	 Determine which line segments are completely inside the disc,
	 and which cross the boundary.
      */
      if(i == seg0) {
	/* initial segment: disc starts from (x0, y0) */
	allin = covered[A] && covered[B];
	bdry  = !allin;
	if(bdry) {
	  if(!covered[A]) nends++;
	  if(!covered[B]) nends++;
	}
      } else {
	/* another segment: disc extends in from either endpoint */
	fromi = from[i];
	toi   = to[i];
	reachable = (covered[fromi] || covered[toi]);
	if(reachable) {
	  allin = (resid[fromi] + resid[toi] >= lengths[i]);
	  bdry = !allin;
	} else allin = bdry = FALSE;
	if(bdry) {
	  if(covered[fromi]) nends++;
	  if(covered[toi]) nends++;
	}
      }
    }
  nendpoints[m] = nends;
  }
}
