/* 
   Definitions of types and data structures for Metropolis-Hastings

   State       Current state of point pattern

   Model       Model parameters passed from R

   Algor       Algorithm parameters (p, q, nrep etc)

   Propo       Proposal in Metropolis-Hastings algorithm

   Cifns       Set of functions for computing the conditional intensity
               for a point process model. 
	       This consists of three functions
                    init(State, Model, Algor) .... initialises auxiliary data
		    eval(State, Propo) ........... evaluates cif
		    update(State,Propo) .......... updates auxiliary data
 */

/* Current state of point pattern */
typedef struct State { 
  double *x;     /* vectors of Cartesian coordinates */
  double *y;
  int *marks;    /* vector of mark values */
  int npts;       /* current number of points */
  int npmax;      /* storage limit */
} State;

/* Parameters of model passed from R */
typedef struct Model {
  double *par;     /* vector of model parameters */
  double *period;  /* width & height of rectangle, if torus */
  int ntypes;      /* number of possible marks */
} Model;

/* RMH Algorithm parameters */
typedef struct Algor {
  double p;         /* probability of proposing shift */
  double q;         /* conditional probability of proposing death */
  int fixall;       /* if TRUE, only shifts of location are feasible */
  int ncond;        /* For conditional simulation, 
		       the first 'ncond' points are fixed */
  int nrep;        /* number of iterations */
  int nverb;        /* print report every 'nverb' iterations */
} Algor;

/* Metropolis-Hastings proposal */
typedef struct Propo {
  double u;         /* location of point of interest */
  double v;
  int mrk;       /* mark of point of interest */
  int ix;           /* index of point of interest, if already in pattern */
  int itype;        /* transition type */
} Propo;

/* transition codes 'itype' */
#define REJECT 0
#define BIRTH 1
#define DEATH 2
#define SHIFT 3

/* conditional intensity functions */

typedef void   (*initfunptr)(State state, Model model, Algor algo);
typedef double (*evalfunptr)(Propo prop, State state);
typedef void   (*updafunptr)(State state, Propo prop);

typedef struct Cifns {
  initfunptr init;
  evalfunptr eval;
  updafunptr update;
  int        marked;
} Cifns;

#define NEED_UPDATE(X) ((X).update != (updafunptr) NULL)

#define NULL_CIFNS { (initfunptr) NULL, (evalfunptr) NULL, (updafunptr) NULL, FALSE}

/* miscellaneous macros */

#ifndef TRUE
#define TRUE (0==0)
#endif
#ifndef FALSE
#define FALSE (!TRUE)
#endif

# define MAT(X,I,J,M) (X[(I)+(J)*(M)])




