/*

   done.c

   $Revision: 1.4 $   $Date: 2006/08/22 02:23:09 $

   Code by Dominic Schuhmacher

   Modified by Adrian Baddeley

*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>

typedef struct State {
  int n; 
  /* vectors of length n */
  int *assig;            /* assignment */
  int *rowlab, *collab;  /* row and col labels */
  int *helper;           /* helping vector to store intermediate results */
  /* n by n matrices */
  int *d;                /* matrix of costs */              
  int *collectvals;    
} State;

#define COST(I,J,STATE,NVALUE) ((STATE)->d)[(NVALUE) * (J) + (I)]

int arraymin(int *a, int n);
void initcosts(State *state);
void maxflow(State *state);
void updatecosts(State *state);
unsigned char colnumassigned(int j, State *state);
void reassign(int startcol, State *state);

/* ------------ The main function ----------------------------- */

void done_R(int *d, int *num, int *assignment)
{
   int i,j; /* indices */
   int n;
   unsigned char feasible = 0; /* boolean for main loop */
   State state;

   /* inputs */
   state.n = n = *num;
   state.d = d;
   /* scratch space */
   state.assig = (int *) R_alloc((long) n, sizeof(int));
   state.rowlab = (int *) R_alloc((long) n, sizeof(int));
   state.collab = (int *) R_alloc((long) n, sizeof(int));
   state.helper = (int *) R_alloc((long) n, sizeof(int));
   state.collectvals = (int *) R_alloc((long) (n * n), sizeof(int));

/* Initialize costm */
   initcosts(&state);

/* For testing: print out cost matrix 
   for (i = 0; i < n; ++i) {
   for (j = 0; j < n; ++j) {
      printf("%d ", COSTM(i, j, &state, n);
   }
   printf("\n");
   }                                  */

/* The main loop */	 
   while(feasible == 0) {
      maxflow(&state);
      if (arraymin(state.assig, n) == -1) {
         updatecosts(&state);
      }
      else {
         feasible = 1;
      }
   }

/* "Return" the final assignment */
   for (i = 0; i < n; i++) {
      assignment[i] = state.assig[i] + 1;
   }

}




/* ------------ Functions called by done_R ------------------------- */


/* Minimal element of an integer array */
int arraymin(int *a, int n) {
  int i, amin;
  if(n < 1)
    return(-1);
  amin = a[0];
  if(n > 1)
    for(i = 0; i < n; i++)
      if(a[i] < amin) amin = a[i];
  return(amin);
}


/* Initialize cost matrix: subtract in each row its minimal entry (from all the
entries in the row), then subtract in each column its minimal entry (from all the
entries in the column) */
void initcosts(State *state) {
   int i,j,m,n;

   n = state->n;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) 
         state->helper[j] = COST(i, j, state, n);
      m = arraymin(state->helper, n);
      for (j = 0; j < n; j++) 
         COST(i, j, state, n) -= m;
   }
   for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) 
	state->helper[i] = COST(i, j, state, n);
      m = arraymin(state->helper, n);
      for (i = 0; i < n; i++) 
	COST(i, j, state, n) -= m;
   }
}


/* Maximize the flow on the (zeros of the) current cost matrix */
void maxflow(State *state) {
   int breakthrough; /* col. no. in which breakthrough occurs */
   unsigned char labelfound = 1; /* 0 if no more labels can be found */
   int i,j,n;

   n = state->n;
   
/* determination of initial assignment for updated cost matrix */
/* Note: this step is not necessary (not used in java version) provided
/* we initialize the assig[i]'s by -1 at the very beginning! */
/* but it seems that it makes the program faster (at least in R) */
   for (i = 0; i < n; i++) 
     state->assig[i] = -1;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (COST(i,j,state,n) == 0 && colnumassigned(j, state) == 0) {
            state->assig[i] = j;
            break;
         }
      }
   }
   while (labelfound == 1) {
      breakthrough = -1;
      /* initialize labels */
      for (i = 0; i < n; i++) {
         state->rowlab[i] = (state->assig[i] == -1) ? -5 : -1;
         state->collab[i] = -1;
         /* -1 means "no index", -5 means "source label" (rows only) */
         /* Note: I think it is not necessary to "properly define" the assignment
            in the first round: everywhere -1 is enough  */
      }
      while (labelfound == 1 && breakthrough == -1) {
         labelfound = 0;
         /* label unlabeled column j that permits flow from some labeled row i */
         /* ("permits flow" means costm[i][j] = 0). Do so for every j          */
         for (i = 0; i < n; i++) {
            if (state->rowlab[i] != -1) {
               for (j = 0; j < n; j++) {
                  if (COST(i, j, state, n) == 0 && state->collab[j] == -1) {
                     state->collab[j] = i;
                     labelfound = 1;
                     if (colnumassigned(j, state) == 0 && breakthrough == -1)
                        breakthrough = j;
                  }
               }
            }
         }
         /* label unlabeled row i that already sends flow to some labeled col j */
         /* ("already sends" means assig[i] = j). Do so for every i             */
         for (j = 0; j < n; j++) {
            if (state->collab[j] != -1) {
               for (i = 0; i < n; i++) {
                  if (state->assig[i] == j && state->rowlab[i] == -1) {
                     state->rowlab[i] = j;
                     labelfound = 1;
                  }
               }
            }
         }
      }
      if (breakthrough != -1) reassign(breakthrough, state);
   }
}


/* Update the cost matrix (called if solution of restricted primal is not feasible
for the original problem): determine the minimum over the submatrix given by all
labeled rows and unlabeled columns, and subtract it from all labeled rows and add
it to all labeled columns. */
void updatecosts(State *state) 
{
   int i,j, n, mini;
   int count = 0; 

   n = state->n;

   for (i = 0; i < n; i++) {
     for (j = 0; j < n; j++) {
       if (state->rowlab[i] != -1 && state->collab[j] == -1) {
	 state->collectvals[count] = COST(i, j, state, n);
	 count++;
       }
     }
   }
   mini = arraymin(state->collectvals, count);
   for (i = 0; i < n; i++) {
     if (state->rowlab[i] != -1) {
       for (j = 0; j < n; j++)
	 COST(i, j, state, n) -= mini;
     }
   }
   for (j = 0; j < n; j++){
     if (state->collab[j] != -1) {
       for (i = 0; i < n; i++) {
	 COST(i, j, state, n) += mini;
       }
     }
   }
}


/* Return 1 if column j is already assigned and 0 otherwise */
unsigned char colnumassigned(int j, State *state) {
   int i, n;

   n = state->n;

   for (i = 0; i < n; i++) {
      if (state->assig[i] == j) return(1);
   }
   return(0);
}


/* Reassign points of one point pattern and points of the other point pattern
according to the row and column labels starting in column startcol */
void reassign(int startcol, State *state) {
   int k,l;

   l = startcol;
   while (l != -5) {
      k = state->collab[l];
      state->assig[k] = l;
      l = state->rowlab[k];
   }
}
