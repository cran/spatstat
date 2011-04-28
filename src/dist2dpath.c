#include <R.h>

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


void dist2dpath(nv, d, adj, dpath, niter) 
  int *nv;     /* number of vertices */
  double *d;  /* matrix of edge lengths */
  int *adj;   /* 0/1 edge matrix of graph */
  double *dpath; /* output - shortest path distance matrix */
  int *niter;
{
  int i, j, k, n, iter;
  int changed, changedij;
  double dij, dik, dkj, dikj;

  n = *nv;

  /* initialise */
  for(i = 0; i < n; i++) 
    for(j = 0; j < n; j++) 
      DPATH(i, j) = (i == j) ? 0 : ((ADJ(i,j)) ? D(i, j) : INFIN);

  /* */
  for(iter = 0; iter < n; iter++) {
    changed = 0;
#ifdef DEBUG
    Rprintf("--------- iteration %d ---------------\n", iter);
#endif
    for(i = 1; i < n; i++) {
      for(j = 0; j < i; j++) {
	dij = DPATH(i,j);
#ifdef DEBUG
	Rprintf("i=%d j=%d dij=%lf\n", i, j, dij);
#endif
	changedij = 0;
	for(k = 0; k < n; k++) {
	  dik = DPATH(i,k);
	  dkj = DPATH(k,j);
	  if(FINITE(dik) && FINITE(dkj)) {
	    dikj = dik + dkj;
#ifdef DEBUG
	    Rprintf("considering %d -> %d -> %d, \t dikj=%lf\n", 
		    i, j, k, dikj);
#endif
	    if(!FINITE(dij) || dikj < dij) {
	      changedij = 1;
#ifdef DEBUG
	      Rprintf("updating i=%d j=%d via k=%d from %lf to %lf\n", 
		      i, j, k, dij, dikj);
#endif
	      dij = dikj;
	    }
	  }
	}
	if(changedij != 0) {
#ifdef DEBUG
	  Rprintf("Resetting d(%d,%d) to %lf\n", i, j, dij);
#endif
	  DPATH(i,j) = DPATH(j,i) = dij;
	  changed = 1;
	}
      }
    }
    if(changed == 0) 
      break;
  }
  
  *niter = iter;
}

