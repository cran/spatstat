#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "methas.h"
#include "chunkloop.h"

void fexitc(const char *msg);


/* To switch on debugging code, 
   insert the line: #define MH_DEBUG TRUE
*/
#ifndef MH_DEBUG
#define MH_DEBUG FALSE
#endif

/* 
   This is the value of 'ix' when we are proposing a birth.
   It must be equal to -1 so that NONE+1 = 0.
*/
#define NONE -1

extern Cifns getcif(char *);

SEXP xmethas(
	     SEXP ncif,
	     SEXP cifname,
	     SEXP beta,
	     SEXP ipar,
	     SEXP iparlen,
	     SEXP period,
	     SEXP xprop,
	     SEXP yprop,
	     SEXP mprop,
	     SEXP ntypes,
	     SEXP nrep,
	     SEXP p,
	     SEXP q,
	     SEXP nverb,
	     SEXP nrep0,
	     SEXP x,
	     SEXP y,
	     SEXP marks,
	     SEXP ncond,
	     SEXP fixall,
             SEXP track)
{
  char *cifstring;
  double cvd, cvn, qnodds, anumer, adenom, betavalue;
  double *iparvector;
  int verb, marked, mustupdate, itype;
  int nfree;
  int irep, ix, j, maxchunk, iverb;
  int Ncif; 
  int *plength;
  long Nmore;
  double *xx, *yy, *xpropose, *ypropose;
  int    *mm,      *mpropose, *pp, *aa;
  SEXP out, xout, yout, mout, pout, aout;
  int tracking;

  State state;
  Model model;
  Algor algo;
  Propo birthprop, deathprop, shiftprop;
  History history;

  /* The following variables are used only for a non-hybrid interaction */
  Cifns thecif;     /* cif structure */
  Cdata *thecdata;  /* pointer to initialised cif data block */

  /* The following variables are used only for a hybrid interaction */
  Cifns *cif;       /* vector of cif structures */
  Cdata **cdata;    /* vector of pointers to initialised cif data blocks */
  int *needupd;     /* vector of logical values */
  int   k;          /* loop index for cif's */

  /* =================== Protect R objects from garbage collector ======= */

  PROTECT(ncif      = AS_INTEGER(ncif)); 
  PROTECT(cifname   = AS_CHARACTER(cifname)); 
  PROTECT(beta      = AS_NUMERIC(beta)); 
  PROTECT(ipar      = AS_NUMERIC(ipar)); 
  PROTECT(iparlen   = AS_INTEGER(iparlen)); 
  PROTECT(period    = AS_NUMERIC(period)); 
  PROTECT(xprop     = AS_NUMERIC(xprop)); 
  PROTECT(yprop     = AS_NUMERIC(yprop)); 
  PROTECT(mprop     = AS_INTEGER(mprop)); 
  PROTECT(ntypes    = AS_INTEGER(ntypes)); 
  PROTECT(nrep      = AS_INTEGER(nrep)); 
  PROTECT(   p      = AS_NUMERIC(p)); 
  PROTECT(   q      = AS_NUMERIC(q)); 
  PROTECT(nverb     = AS_INTEGER(nverb)); 
  PROTECT(nrep0     = AS_INTEGER(nrep0)); 
  PROTECT(   x      = AS_NUMERIC(x)); 
  PROTECT(   y      = AS_NUMERIC(y)); 
  PROTECT( marks    = AS_INTEGER(marks)); 
  PROTECT(fixall    = AS_INTEGER(fixall)); 
  PROTECT(ncond     = AS_INTEGER(ncond)); 
  PROTECT(track     = AS_INTEGER(track)); 

                    /* that's 21 protected objects */

  /* =================== Translate arguments from R to C ================ */

  /* 
     Ncif is the number of cif's
     plength[i] is the number of interaction parameters in the i-th cif
  */
  Ncif = *(INTEGER_POINTER(ncif));
  plength = INTEGER_POINTER(iparlen);

  /* copy RMH algorithm parameters */
  algo.nrep   = *(INTEGER_POINTER(nrep));
  algo.nverb  = *(INTEGER_POINTER(nverb));
  algo.nrep0  = *(INTEGER_POINTER(nrep0));
  algo.p = *(NUMERIC_POINTER(p));
  algo.q = *(NUMERIC_POINTER(q));
  algo.fixall = ((*(INTEGER_POINTER(fixall))) == 1);
  algo.ncond =  *(INTEGER_POINTER(ncond));

  /* copy model parameters without interpreting them */
  model.beta = NUMERIC_POINTER(beta);
  model.ipar = iparvector = NUMERIC_POINTER(ipar);
  model.period = NUMERIC_POINTER(period);
  model.ntypes = *(INTEGER_POINTER(ntypes));
  marked = (model.ntypes > 1);
  
  /* copy initial state */
  state.npts   = LENGTH(x);
  state.npmax  = 4 * ((state.npts > 256) ? state.npts : 256);
  state.x = (double *) R_alloc(state.npmax, sizeof(double));
  state.y = (double *) R_alloc(state.npmax, sizeof(double));
  xx = NUMERIC_POINTER(x);
  yy = NUMERIC_POINTER(y);
  if(marked) {
    state.marks =(int *) R_alloc(state.npmax, sizeof(int));
    mm = INTEGER_POINTER(marks);
  }
  if(!marked) {
    for(j = 0; j < state.npts; j++) {
      state.x[j] = xx[j];
      state.y[j] = yy[j];
    }
  } else {
    for(j = 0; j < state.npts; j++) {
      state.x[j] = xx[j];
      state.y[j] = yy[j];
      state.marks[j] = mm[j];
    }
  }
#if MH_DEBUG
  Rprintf("\tnpts=%d\n", state.npts);
#endif

  /* access proposal data */
  xpropose = NUMERIC_POINTER(xprop);
  ypropose = NUMERIC_POINTER(yprop);
  mpropose = INTEGER_POINTER(mprop);
  /* we need to initialise 'mpropose' to keep compilers happy.
     mpropose is only used for marked patterns.
     Note 'mprop' is always a valid pointer */

  
  /* ================= Allocate space for cifs etc ========== */

  if(Ncif > 1) {
    cif = (Cifns *) R_alloc(Ncif, sizeof(Cifns));
    cdata = (Cdata **) R_alloc(Ncif, sizeof(Cdata *));
    needupd = (int *) R_alloc(Ncif, sizeof(int));
  } else {
    /* Keep the compiler happy */
    cif = (Cifns *) R_alloc(1, sizeof(Cifns));
    cdata = (Cdata **) R_alloc(1, sizeof(Cdata *));
    needupd = (int *) R_alloc(1, sizeof(int));
  }


  /* ================= Determine process to be simulated  ========== */
  
  /* Get the cif's */
  if(Ncif == 1) {
    cifstring = (char *) STRING_VALUE(cifname);
    thecif = getcif(cifstring);
    mustupdate = NEED_UPDATE(thecif);
    if(thecif.marked && !marked)
      fexitc("cif is for a marked point process, but proposal data are not marked points; bailing out.");
    /* Keep compiler happy*/
    cif[0] = thecif;
    needupd[0] = mustupdate;
  } else {
    mustupdate = FALSE;
    for(k = 0; k < Ncif; k++) {
      cifstring = (char *) CHAR(STRING_ELT(cifname, k));
      cif[k] = getcif(cifstring);
      needupd[k] = NEED_UPDATE(cif[k]);
      if(needupd[k])
	mustupdate = TRUE;
      if(cif[k].marked && !marked)
	fexitc("component cif is for a marked point process, but proposal data are not marked points; bailing out.");
    }
  }
  /* ============= Initialise transition history ========== */

  tracking = (*(INTEGER_POINTER(track)) != 0);
  if(tracking) {
    history.nmax = algo.nrep;
    history.n = 0;
    history.proptype = (int *) R_alloc(algo.nrep, sizeof(int));
    history.accepted = (int *) R_alloc(algo.nrep, sizeof(int));
  }

  /* ================= Initialise algorithm ==================== */
 
  /* Interpret the model parameters and initialise auxiliary data */
  if(Ncif == 1) {
    thecdata = (*(thecif.init))(state, model, algo);
    /* keep compiler happy */
    cdata[0] = thecdata;
  } else {
    for(k = 0; k < Ncif; k++) {
      if(k > 0)
	model.ipar += plength[k-1];
      cdata[k] = (*(cif[k].init))(state, model, algo);
    }
  }

  /* Set the fixed elements of the proposal objects */
  birthprop.itype = BIRTH;
  deathprop.itype = DEATH;
  shiftprop.itype = SHIFT;
  birthprop.ix = NONE;
  if(!marked) 
    birthprop.mrk = deathprop.mrk = shiftprop.mrk = NONE;

  /* Set up some constants */
  verb   = (algo.nverb !=0);
  qnodds = (1.0 - algo.q)/algo.q;


  /* Set value of beta for unmarked process */
  /* (Overwritten for marked process, but keeps compiler happy) */
  betavalue = model.beta[0];

  /* ============= Run Metropolis-Hastings  ================== */

  /* Initialise random number generator */
  GetRNGstate();

  if(tracking) {
    /* saving transition history */
#define MH_TRACKING TRUE
    if(Ncif == 1) {
      /* single interaction */
#define MH_SINGLE TRUE
      if(marked) {
	/* marked process */
#define MH_MARKED TRUE

	/* run loop */
#include "mhloop.h"

      } else {
	/* unmarked process */
#undef MH_MARKED
#define MH_MARKED FALSE

	/* run loop */
#include "mhloop.h"

      }
    } else {
      /* hybrid interaction */
#undef MH_SINGLE
#define MH_SINGLE FALSE
      if(marked) {
	/* marked process */
#undef MH_MARKED
#define MH_MARKED TRUE

	/* run loop */
#include "mhloop.h"

      } else {
	/* unmarked process */
#undef MH_MARKED
#define MH_MARKED FALSE

	/* run loop */
#include "mhloop.h"

      }
    }
  } else {
    /* not saving transition history */
#undef MH_TRACKING
#define MH_TRACKING FALSE
    if(Ncif == 1) {
      /* single interaction */
#undef MH_SINGLE
#define MH_SINGLE TRUE
      if(marked) {
	/* marked process */
#undef MH_MARKED
#define MH_MARKED TRUE

	/* run loop */
#include "mhloop.h"

      } else {
	/* unmarked process */
#undef MH_MARKED
#define MH_MARKED FALSE

	/* run loop */
#include "mhloop.h"

      }
    } else {
      /* hybrid interaction */
#undef MH_SINGLE
#define MH_SINGLE FALSE
      if(marked) {
	/* marked process */
#undef MH_MARKED
#define MH_MARKED TRUE

	/* run loop */
#include "mhloop.h"

      } else {
	/* unmarked process */
#undef MH_MARKED
#define MH_MARKED FALSE

	/* run loop */
#include "mhloop.h"

      }
    }
  }

  /* relinquish random number generator */
  PutRNGstate();

  /* ============= Done  ================== */

  /* Create space for output, and copy final state */
  /* Point coordinates */
  PROTECT(xout = NEW_NUMERIC(state.npts));
  PROTECT(yout = NEW_NUMERIC(state.npts));
  xx = NUMERIC_POINTER(xout);
  yy = NUMERIC_POINTER(yout);
  for(j = 0; j < state.npts; j++) {
    xx[j] = state.x[j];
    yy[j] = state.y[j];
  }
  /* Marks */
  if(marked) {
    PROTECT(mout = NEW_INTEGER(state.npts));
    mm = INTEGER_POINTER(mout);
    for(j = 0; j < state.npts; j++) 
      mm[j] = state.marks[j];
  } else {
    /* Keep the compiler happy */
    PROTECT(mout = NEW_INTEGER(1));
    mm = INTEGER_POINTER(mout);
    mm[0] = 0;
  }
  /* Transition history */
  if(tracking) {
    PROTECT(pout = NEW_INTEGER(algo.nrep));
    PROTECT(aout = NEW_INTEGER(algo.nrep));
    pp = INTEGER_POINTER(pout);
    aa = INTEGER_POINTER(aout);
    for(j = 0; j < algo.nrep; j++) {
      pp[j] = history.proptype[j];
      aa[j] = history.accepted[j];
    }
  } else {
    /* Keep the compiler happy */
    PROTECT(pout = NEW_INTEGER(1));
    PROTECT(aout = NEW_INTEGER(1));
    pp = INTEGER_POINTER(pout);
    aa = INTEGER_POINTER(aout);
    pp[0] = aa[0] = 0;
  }

  /* Pack up into list object for return */
  if(!tracking) {
    /* no transition history */
    if(!marked) {
      PROTECT(out = NEW_LIST(2));
      SET_VECTOR_ELT(out, 0, xout);
      SET_VECTOR_ELT(out, 1, yout);
    } else {
      PROTECT(out = NEW_LIST(3)); 
      SET_VECTOR_ELT(out, 0, xout);
      SET_VECTOR_ELT(out, 1, yout); 
      SET_VECTOR_ELT(out, 2, mout);
    }
  } else {
    /* transition history */
    if(!marked) {
      PROTECT(out = NEW_LIST(4));
      SET_VECTOR_ELT(out, 0, xout);
      SET_VECTOR_ELT(out, 1, yout);
      SET_VECTOR_ELT(out, 2, pout);
      SET_VECTOR_ELT(out, 3, aout);
      } else {
      PROTECT(out = NEW_LIST(5)); 
      SET_VECTOR_ELT(out, 0, xout);
      SET_VECTOR_ELT(out, 1, yout); 
      SET_VECTOR_ELT(out, 2, mout);
      SET_VECTOR_ELT(out, 3, pout);
      SET_VECTOR_ELT(out, 4, aout);
    }
  }
  UNPROTECT(27);  /* 21 arguments plus xout, yout, mout, pout, aout, out */
  return(out);
}
