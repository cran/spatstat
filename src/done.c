/*

   done.c

   $Revision: 1.1 $   $Date: 2006/06/30 09:53:52 $

   Code by Dominic Schuhmacher
  
*/

#include <stdio.h>
#include <stdlib.h>

int intmin(int a, int b);
int arraymin(int *a, int n);
void initcosts();
void maxflow();
void updatecosts();
unsigned char colnumassigned(int j);
void reassign(int startcol);

/* external variables */
int n;
int *assig, *rowlab, *collab, *helper; /* assignment, row and col labels */
                 /* helper = helping vector to store intermediate results */
int **costm;

/* ----------------------------- */
/* --- Version of 21/02/2006 --- */
/* ----------------------------- */

/* ------------ The main function ----------------------------- */

void done_R(int *d, int *num, int *assignment) {
   int i,j; /* indices */
   unsigned char feasible = 0; /* boolean for main loop */

   n = *num;

/* dynamic allocation of rowlab, collab, costm */
   helper = (int *)malloc(n * sizeof(int));
   rowlab = (int *)malloc(n * sizeof(int));
   collab = (int *)malloc(n * sizeof(int));
   assig = (int *)malloc(n * sizeof(int));
   if (helper == NULL || rowlab == NULL || collab == NULL || assig == NULL) {
      printf("out of memory\n");
      exit(1);
   }
   costm = malloc(n * sizeof(int *));
   if (costm == NULL) {
      printf("out of memory\n");
      exit(1);
   }
   for (i = 0; i < n; i++) {
      costm[i] = malloc(n * sizeof(int));
      if (costm[i] == NULL) {
         printf("out of memory\n");
         exit(1);
      }
   }

/* Copy distance matrix obtained from R into costm */
   for (i = 0; i < n; ++i) {
   for (j = 0; j < n; ++j) {
      costm[i][j] = d[n * j + i];
   }
   }

/* Initialize costm */
   initcosts();

/* For testing: print out cost matrix 
   for (i = 0; i < n; ++i) {
   for (j = 0; j < n; ++j) {
      printf("%d ", costm[i][j]);
   }
   printf("\n");
   }                                  */

/* The main loop */	 
   while(feasible == 0) {
      maxflow();
      if (arraymin(assig, n) == -1) {
         updatecosts();
      }
      else {
         feasible = 1;
      }
   }

/* "Return" the final assignment */
   for (i = 0; i < n; i++) {
      assignment[i] = assig[i] + 1;
   }

/* Free the memory used by the dyn. allocated stuff */
   for(i = 0; i < n; i++) {
      free((void *)costm[i]);
   }
   free((void *)costm);
   free(assig);
   free(collab);
   free(rowlab);
   free(helper);
}




/* ------------ Functions called by done_R ------------------------- */


/* Minimum of two integers */
int intmin(int a, int b) {
   return(a < b ? a : b);
}


/* Minimal element of an array */
int arraymin(int *a, int n) {
   if (n == 1) return(a[0]);
   else return(intmin(arraymin(a, n-1), a[n-1]));
}


/* Initialize cost matrix: subtract in each row its minimal entry (from all the
entries in the row), then subtract in each column its minimal entry (from all the
entries in the column) */
void initcosts() {
   int i,j;
   int m;

   for (i = 0; i < n; i++) {
      m = arraymin(costm[i], n);
      for (j = 0; j < n; j++) {
         costm[i][j] = costm[i][j] - m;
      }
   }
   for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
         helper[i] = costm[i][j];
      }
      m = arraymin(helper, n);
      for (i = 0; i < n; i++) {
         costm[i][j] = costm[i][j] - m;
      }
   }
}


/* Maximize the flow on the (zeros of the) current cost matrix */
void maxflow() {
   int breakthrough; /* col. no. in which breakthrough occurs */
   unsigned char labelfound = 1; /* 0 if no more labels can be found */
   int i,j;

/* determination of initial assignment for updated cost matrix */
/* Note: this step is not necessary (not used in java version) provided
/* we initialize the assig[i]'s by -1 at the very beginning! */
/* but it seems that it makes the program faster (at least in R) */
   for (i = 0; i < n; i++) assig[i] = -1;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (costm[i][j] == 0 && colnumassigned(j) == 0) {
            assig[i] = j;
            break;
         }
      }
   }
   while (labelfound == 1) {
      breakthrough = -1;
      /* initialize labels */
      for (i = 0; i < n; i++) {
         rowlab[i] = (assig[i] == -1) ? -5 : -1;
         collab[i] = -1;
         /* -1 means "no index", -5 means "source label" (rows only) */
         /* Note: I think it is not necessary to "properly define" the assignment
            in the first round: everywhere -1 is enough  */
      }
      while (labelfound == 1 && breakthrough == -1) {
         labelfound = 0;
         /* label unlabeled column j that permits flow from some labeled row i */
         /* ("permits flow" means costm[i][j] = 0). Do so for every j          */
         for (i = 0; i < n; i++) {
            if (rowlab[i] != -1) {
               for (j = 0; j < n; j++) {
                  if (costm[i][j] == 0 && collab[j] == -1) {
                     collab[j] = i;
                     labelfound = 1;
                     if (colnumassigned(j) == 0 && breakthrough == -1)
                        breakthrough = j;
                  }
               }
            }
         }
         /* label unlabeled row i that already sends flow to some labeled col j */
         /* ("already sends" means assig[i] = j). Do so for every i             */
         for (j = 0; j < n; j++) {
            if (collab[j] != -1) {
               for (i = 0; i < n; i++) {
                  if (assig[i] == j && rowlab[i] == -1) {
                     rowlab[i] = j;
                     labelfound = 1;
                  }
               }
            }
         }
      }
      if (breakthrough != -1) reassign(breakthrough);
   }
}


/* Update the cost matrix (called if solution of restricted primal is not feasible
for the original problem): determine the minimum over the submatrix given by all
labeled rows and unlabeled columns, and subtract it from all labeled rows and add
it to all labeled columns. */
void updatecosts() {
   int i,j, mini;
   int count = 0; 
   int *collectvals; /* probably program is a little faster if this is defined externally */

   collectvals = (int *)malloc(n * n * sizeof(int));
   if (collectvals == NULL) {
      printf("out of memory\n");
      exit(1);
   }

   for (i = 0; i < n; i++) {
   for (j = 0; j < n; j++) {
      if (rowlab[i] != -1 && collab[j] == -1) {
         collectvals[count] = costm[i][j];
         count++;
      }
   }
   }
   mini = arraymin(collectvals, count);
   for (i = 0; i < n; i++) {
      if (rowlab[i] != -1) {
         for (j = 0; j < n; j++) {
            costm[i][j] -= mini;
         }
      }
   }
   for (j = 0; j < n; j++){
      if (collab[j] != -1) {
         for (i = 0; i < n; i++) {
            costm[i][j] += mini;
         }
      }
   }
   free(collectvals);
}


/* Return 1 if column j is already assigned and 0 otherwise */
unsigned char colnumassigned(int j) {
   int i;

   for (i = 0; i < n; i++) {
      if (assig[i] == j) return(1);
   }
   return(0);
}


/* Reassign points of one point pattern and points of the other point pattern
according to the row and column labels starting in column startcol */
void reassign(int startcol) {
   int k,l;

   l = startcol;
   while (l != -5) {
      k = collab[l];
      assig[k] = l;
      l = rowlab[k];
   }
}
