/*

   dinfty.c

   $Revision: 1.1 $   $Date: 2006/06/30 09:50:00 $

   Code by Dominic Schuhmacher

*/

#include <stdio.h>
#include <stdlib.h>

int intmax(int a, int b);
int arraymax(int *a, int n);
void swap(int i, int j);
int largestmobpos();

/* external variables */
int n, currmin;
int *current, *travel, *mobile; 
int **costm;

/* ----------------------------- */
/* --- Version of 01/03/2006 --- */
/* ----------------------------- */

/* ------------ The main function ----------------------------- */

void dinfty_R(int *d, int *num, int *assignment) {
   int *assig, *distrelev; 
   int i,j; /* indices */
   int lmp, lmq; /* largest mobile position and its neighbor */
   int oldcost, newmax;

   n = *num;

/* dynamic allocation of travel, mobile, current, and assig */
   travel = (int *)malloc(n * sizeof(int));
   mobile = (int *)malloc(n * sizeof(int));
   current = (int *)malloc(n * sizeof(int));
   assig = (int *)malloc(n * sizeof(int));
   distrelev = (int *)malloc(n * sizeof(int));
   if (travel == NULL || mobile == NULL || current == NULL || assig == NULL || distrelev == NULL) {
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
      if(costm[i] == NULL) {
         printf("out of memory\n");
         exit(1);
      }
   }

/* Copy distance matrix obtained from R into costm */
   for (i = 0; i < n; i++) {
   for (j = 0; j < n; j++) {
      costm[i][j] = d[n * j + i];
   }
   }

/*                                                               */
/* We use the Johnson-Trotter Algorithm for listing permutations */
/*                                                               */

/* Initialize the algorithm */
   for (i = 0; i < n; i++) {
      travel[i] = -1;   /* all numbers traveling to the left */
      mobile[i] = 1;    /* all numbers mobile */
      current[i] = i;   /* current permutation is the identity */
      assig[i] = i;     /* best permutation up to now is the identity */
      distrelev[i] = costm[i][i];   /* pick relevant entries in the cost matrix */
   }
   currmin = arraymax(distrelev, n);   /* minimal max up to now */

/* The main loop */
   while(arraymax(mobile, n) == 1) {
      lmp = largestmobpos();
      lmq = lmp + travel[lmp];
      swap(lmp, lmq);
      for (i = 0; i < n; i++) {
         if (current[i] > current[lmq])
            travel[i] = -travel[i];
         j = i + travel[i];
         if (j < 0 || j > n-1 || current[i] < current[j])
            mobile[i] = 0;
         else
            mobile[i] = 1;
         distrelev[i] = costm[i][current[i]];
      }
      /* Calculation of new maximal value */
      newmax = arraymax(distrelev, n);
      if (newmax < currmin) {
         currmin = newmax;
         for (i = 0; i < n; i++) {
            assig[i] = current[i];
         }
      }
   }
/* For testing: print distance from within C program
   printf("Prohorov distance is %d\n", currmin);     /*

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
   free(distrelev);
   free(current);
   free(travel);
   free(mobile);
}


/* Maximum of two integers */
int intmax(int a, int b) {
   return(a > b ? a : b);
}


/* Maximal element of an array */
int arraymax(int *a, int n) {
   if (n == 1) return(a[0]);
   else return(intmax(arraymax(a, n-1), a[n-1]));
}


/* Swap elements i and j in current and in travel */
void swap(int i, int j) {
   int v;

   v = current[i];
   current[i] = current[j];
   current[j] = v;
   v = travel[i];
   travel[i] = travel[j];
   travel[j] = v;
}


/* Return index of largest mobile number in current */
int largestmobpos() {
   int i,j, maxval;
   int *collectvals; /* probably program is a little faster if this is defined externally */

   collectvals = (int *)malloc(n * n * sizeof(int));
   if (collectvals == NULL) {
      printf("out of memory\n");
      exit(1);
   }
   j = 0;
   for (i = 0; i < n; i++) {
      if (mobile[i] == 1) {
         collectvals[j] = current[i];
         j++;
      }
   }
   maxval = arraymax(collectvals, j);
   for (i = 0; i < n; i++) {
      if (current[i] == maxval) {
         free(collectvals);
         return(i);
      }
   }
   printf("something's fishy...\n");
}
