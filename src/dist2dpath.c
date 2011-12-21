#include <R.h>
#include <R_ext/Utils.h>

/*
  given matrix of edge lengths
  compute matrix of shortest-path distances
*/

#undef DEBUG 

#define MATRIX(X,I,J) (X)[(J) + n * (I)]
#define D(I,J)     MATRIX(d,     I, J)
#define DPATH(I,J) MATRIX(dpath, I, J)
#define ADJ(I,J)   (MATRIX(adj,  I, J) != 0)

#define INFIN -1
#define FINITE(X) ((X) >= 0)

void dist2dpath(nv, d, adj, dpath, tol, niter, status) 
  int *nv;     /* number of vertices */
  double *d;  /* matrix of edge lengths */
  int *adj;   /* 0/1 edge matrix of graph */
  double *tol;  /* tolerance threshold */
  double *dpath; /* output - shortest path distance matrix */
  int *niter, *status; /* status = 0 for convergence */
{
  int i, j, k, n, iter, maxiter, changed;
  double dij, dik, dkj, dikj;
  double eps, diff, maxdiff;

  int totaledges, starti, nneighi, increm, pos;
  int *start, *nneigh, *indx;

  n = *nv;
  eps = *tol;

  /* initialise and count edges */
  *status = -1;
  totaledges = 0;
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      DPATH(i, j) = (i == j) ? 0 : ((ADJ(i,j)) ? D(i, j) : INFIN);
      if((i != j) && ADJ(i,j)) ++totaledges;
    }
  }

  maxiter = 2 + ((totaledges > n) ? totaledges : n);

  /* store indices j for each edge (i,j) */
  indx = (int *) R_alloc(totaledges, sizeof(int));
  nneigh = (int *) R_alloc(n, sizeof(int));
  start  = (int *) R_alloc(n, sizeof(int));

  pos = 0;
  for(i = 0; i < n; i++) {
    nneigh[i] = 0;
    start[i] = pos;
#ifdef DEBUG 
    Rprintf("Neighbours of %d:\n", i);
#endif
    for(j = 0; j < n; j++) {
      if((i != j) && ADJ(i,j) && FINITE(D(i,j))) {
#ifdef DEBUG 
	Rprintf("\t%d\n", j);
#endif
	++(nneigh[i]);
	if(pos > totaledges)
	  error("internal error: pos exceeded storage");
	indx[pos] = j;
	++pos;
      }
    }
  }

  /* run */
  for(iter = 0; iter < maxiter; iter++) {

    changed = 0;
    maxdiff = 0;

#ifdef DEBUG
    Rprintf("--------- iteration %d ---------------\n", iter);
#endif
    for(i = 0; i < n; i++) {
      R_CheckUserInterrupt();
      nneighi = nneigh[i];
      if(nneighi > 0) {
	/* run through neighbours k of i */
	starti = start[i];
	for(increm = 0, pos=starti; increm < nneighi; ++increm, ++pos) {
	  k = indx[pos];
	  dik = DPATH(i,k);
#ifdef DEBUG
	    Rprintf("i=%d k=%d dik=%lf\n", i, k, dik);
#endif
	  /* now run through all other vertices j */
	  for(j = 0; j < n; j++) {
	    if(j != i && j != k) {
	      dij = DPATH(i,j);
	      dkj = DPATH(k,j);
	      if(FINITE(dkj)) {
		dikj = dik + dkj;
#ifdef DEBUG
		Rprintf("considering %d -> (%d) -> %d,\t dij=%lf, dikj=%lf\n", 
			i, k, j, dij, dikj);
#endif
		if(!FINITE(dij) || dikj < dij) {
#ifdef DEBUG
		  Rprintf("updating i=%d j=%d via k=%d from %lf to %lf\n", 
			  i, j, k, dij, dikj);
#endif
		  DPATH(i,j) = DPATH(j,i) = dikj;
		  changed = 1;
		  diff = (FINITE(dij)) ? dij - dikj : dikj;
		  if(diff > maxdiff) maxdiff = diff;
		}
	      }
	    }
	  }
	}
      }
    }
    if(changed == 0) {
      /* algorithm converged */
#ifdef DEBUG
      Rprintf("Algorithm converged\n");
#endif
      *status = 0;
      break;
    } else if(FINITE(maxdiff) && maxdiff < eps) {
      /* tolerance reached */
#ifdef DEBUG
      Rprintf("Algorithm terminated with maxdiff=%lf\n", maxdiff);
#endif
      *status = 1;
      break;
    }
  }

#ifdef DEBUG
  Rprintf("Returning after %d iterations on %d vertices\n", iter, n);
#endif
  
  *niter = iter;
}

