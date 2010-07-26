#include <R.h>
#include <Rdefines.h>
#include "methas.h"

void fexitc(const char *msg);


/* to switch on debugging code, 
   insert the line: #define MHDEBUG 1
*/

/* 
   This is the value of 'ix' when we are proposing a birth.
   It must be equal to -1 so that NONE+1 = 0.
*/
#define NONE -1

extern Cifns getcif(char *);

SEXP xmethas(SEXP cifname,
	    SEXP par,
	    SEXP period,
	    SEXP xprop,
	    SEXP yprop,
	    SEXP mprop,
	    SEXP ntypes,
            SEXP nrep,
	    SEXP p,
	    SEXP q,
	    SEXP nverb,
	    SEXP x,
	    SEXP y,
	    SEXP marks,
	    SEXP ncond,
	    SEXP fixall)
{
  char *cifstring;
  double cvd, cvn, qnodds, anumer, adenom;
  int verb, marked, needupd, itype;
  int nfree;
  int irep, ix, j;
  long Nmore;
  double *xx, *yy, *xpropose, *ypropose;
  int    *mm,      *mpropose;
  SEXP out, xout, yout, mout;

  State state;
  Model model;
  Algor algo;
  Propo birthprop, deathprop, shiftprop;
  Cifns cif;
  Cdata *cdata;

  /* =================== Protect R objects from garbage collector ======= */

  PROTECT(cifname   = AS_CHARACTER(cifname)); 
  PROTECT(par       = AS_NUMERIC(par)); 
  PROTECT(period    = AS_NUMERIC(period)); 
  PROTECT(xprop     = AS_NUMERIC(xprop)); 
  PROTECT(yprop     = AS_NUMERIC(yprop)); 
  PROTECT(mprop     = AS_INTEGER(mprop)); 
  PROTECT(ntypes    = AS_INTEGER(ntypes)); 
  PROTECT(nrep      = AS_INTEGER(nrep)); 
  PROTECT(   p      = AS_NUMERIC(p)); 
  PROTECT(   q      = AS_NUMERIC(q)); 
  PROTECT(nverb     = AS_INTEGER(nverb)); 
  PROTECT(   x      = AS_NUMERIC(x)); 
  PROTECT(   y      = AS_NUMERIC(y)); 
  PROTECT( marks    = AS_INTEGER(marks)); 
  PROTECT(fixall    = AS_INTEGER(fixall)); 
  PROTECT(ncond     = AS_INTEGER(ncond)); 

                    /* that's 16 protected objects */

  /* =================== Translate arguments from R to C ================ */

  cifstring = (char *) STRING_VALUE(cifname);

  /* copy RMH algorithm parameters */
  algo.nrep   = *(INTEGER_POINTER(nrep));
  algo.nverb  = *(INTEGER_POINTER(nverb));
  algo.p = *(NUMERIC_POINTER(p));
  algo.q = *(NUMERIC_POINTER(q));
  algo.fixall = ((*(INTEGER_POINTER(fixall))) == 1);
  algo.ncond =  *(INTEGER_POINTER(ncond));

  /* copy model parameters without interpreting them */
  model.par = NUMERIC_POINTER(par);
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
#ifdef MHDEBUG
  Rprintf("\tnpts=%d\n", state.npts);
#endif

  /* access proposal data */
  xpropose = NUMERIC_POINTER(xprop);
  ypropose = NUMERIC_POINTER(yprop);
  mpropose = INTEGER_POINTER(mprop);
  /* we need to initialise 'mpropose' to keep compilers happy.
     mpropose is only used for marked patterns.
     Note 'mprop' is always a valid pointer */


  /* ================= Determine process to be simulated  ========== */
  
  /* Get the three function pointers */
  cif = getcif(cifstring);

  needupd = NEED_UPDATE(cif);

  if(cif.marked && !marked)
    fexitc("cif is for a marked point process, but proposal data are not marked points; bailing out.");

  /* ================= Initialise algorithm ==================== */
 
  /* Initialise random number generator */
  GetRNGstate();

  /* Interpret the model parameters and initialise auxiliary data */
  cdata = (*(cif.init))(state, model, algo);

  /* set the fixed elements of the proposal objects */
  birthprop.itype = BIRTH;
  deathprop.itype = DEATH;
  shiftprop.itype = SHIFT;
  birthprop.ix = NONE;
  if(!marked) 
    birthprop.mrk = deathprop.mrk = shiftprop.mrk = NONE;

  /* Set up some constants */
  verb   = (algo.nverb !=0);
  qnodds = (1.0 - algo.q)/algo.q;

  /* ============= Run Metropolis-Hastings  ================== */

  for(irep = 0; irep < algo.nrep; irep++) {

#ifdef MHDEBUG
    Rprintf("iteration %d\n", irep);
#endif

    if(verb && ((irep+1)%algo.nverb == 0))
      Rprintf("iteration %d\n", irep+1);

    itype = REJECT;

    nfree = state.npts - algo.ncond;  /* number of 'free' points */

    /* ................  generate proposal ..................... */
    /* Shift or birth/death: */
    if(unif_rand() > algo.p) {
#ifdef MHDEBUG
    Rprintf("propose birth or death\n");
#endif
      /* Birth/death: */
      if(unif_rand() > algo.q) {
	/* Propose birth: */
	birthprop.u = xpropose[irep];
	birthprop.v = ypropose[irep];
	if(marked)
	  birthprop.mrk = mpropose[irep];
#ifdef MHDEBUG
	if(marked)
	  Rprintf("propose birth at (%lf, %lf) with mark %d\n", 
		  birthprop.u, birthprop.v, birthprop.mrk);
	else 
	  Rprintf("propose birth at (%lf, %lf)\n", birthprop.u, birthprop.v);
#endif
	/* evaluate conditional intensity */
	anumer = (*(cif.eval))(birthprop, state, cdata);
	adenom = qnodds*(nfree+1);
#ifdef MHDEBUG
	Rprintf("cif = %lf, Hastings ratio = %lf\n", anumer, anumer/adenom);
#endif
	/* accept/reject */
	if(unif_rand() * adenom < anumer) {
#ifdef MHDEBUG
	  Rprintf("accepted birth\n");
#endif
	  itype = BIRTH;  /* Birth proposal accepted. */
	}
      } else if(nfree > 0) {
	/* Propose death: */
	ix = floor(nfree * unif_rand());
	if(ix < 0) ix = 0;
	ix = algo.ncond + ix;
	if(ix >= state.npts) ix = state.npts - 1;
	deathprop.ix = ix;
	deathprop.u  = state.x[ix];
	deathprop.v  = state.y[ix];
	if(marked) 
	  deathprop.mrk = state.marks[ix];
#ifdef MHDEBUG
	if(marked)
	  Rprintf("propose death of point %d = (%lf, %lf) with mark %d\n", 
		  ix, deathprop.u, deathprop.v, deathprop.mrk);
	else 
	  Rprintf("propose death of point %d = (%lf, %lf)\n", 
		  ix, deathprop.u, deathprop.v);
#endif
	/* evaluate conditional intensity */
	adenom = (*(cif.eval))(deathprop, state, cdata);
	anumer = qnodds * nfree;
#ifdef MHDEBUG
	Rprintf("cif = %lf, Hastings ratio = %lf\n", adenom, anumer/adenom);
#endif
	/* accept/reject */
	if(unif_rand() * adenom < anumer) {
#ifdef MHDEBUG
	  Rprintf("accepted death\n");
#endif
	  itype = DEATH; /* Death proposal accepted. */
	}
      }
    } else if(nfree > 0) {
      /* Propose shift: */
      /* point to be shifted */
      ix = floor(nfree * unif_rand());
      if(ix < 0) ix = 0;
      ix = algo.ncond + ix;
      if(ix >= state.npts) ix = state.npts - 1;
      deathprop.ix = ix;
      deathprop.u  = state.x[ix];
      deathprop.v  = state.y[ix];
      if(marked) 
	deathprop.mrk = state.marks[ix];
      /* where to shift */
      shiftprop.u = xpropose[irep]; 
      shiftprop.v = ypropose[irep];
      if(marked) 
	shiftprop.mrk = (algo.fixall) ? deathprop.mrk : mpropose[irep];
      shiftprop.ix = ix;
#ifdef MHDEBUG
      if(marked)
 	Rprintf("propose shift of point %d = (%lf, %lf)[mark %d] to (%lf, %lf)[mark %d]\n", 
		ix, deathprop.u, deathprop.v, deathprop.mrk, 
		shiftprop.u, shiftprop.v, shiftprop.mrk);
      else 
 	Rprintf("propose shift of point %d = (%lf, %lf) to (%lf, %lf)\n", 
		ix, deathprop.u, deathprop.v, shiftprop.u, shiftprop.v);
#endif
      /* evaluate cif in two stages */
      cvd = (*(cif.eval))(deathprop, state, cdata);
      cvn = (*(cif.eval))(shiftprop, state, cdata);
#ifdef MHDEBUG
	Rprintf("cif[old] = %lf, cif[new] = %lf, Hastings ratio = %lf\n", 
		cvd, cvn, cvn/cvd);
#endif
      /* accept/reject */
      if(unif_rand() * cvd < cvn) {
#ifdef MHDEBUG
	Rprintf("accepted shift\n");
#endif
	itype = SHIFT;          /* Shift proposal accepted . */
      }
    }
    if(itype != REJECT) {
      /* ....... implement the transition ............  */
      if(itype == BIRTH) {      
	/* Birth transition */
	/* add point at (u,v) */
#ifdef MHDEBUG
	if(marked)
	  Rprintf("implementing birth at (%lf, %lf) with mark %d\n", 
		  birthprop.u, birthprop.v, birthprop.mrk);
	else
	  Rprintf("implementing birth at (%lf, %lf)\n", 
		  birthprop.u, birthprop.v);
#endif
	if(state.npts + 1 > state.npmax) {
#ifdef MHDEBUG
	  Rprintf("!!!!!!!!!!! storage overflow !!!!!!!!!!!!!!!!!\n");
#endif
	  /* storage overflow; allocate more storage */
	  Nmore = 2 * state.npmax;
	  state.x = (double *) S_realloc((char *) state.x, 
					 Nmore,  state.npmax, 
					 sizeof(double));
	  state.y = (double *) S_realloc((char *) state.y, 
					 Nmore,  state.npmax, 
					 sizeof(double));
	  if(marked)
	    state.marks = (int *) S_realloc((char *) state.marks, 
					    Nmore,  state.npmax, 
					    sizeof(int));
	  state.npmax = Nmore;

	  /* call the initialiser again, to allocate additional space */
	  cdata = (*(cif.init))(state, model, algo);
	} 
	/* Update auxiliary variables first */
	if(needupd)
	  (*(cif.update))(state, birthprop, cdata);
	/* Now add point */
	state.x[state.npts] = birthprop.u;
	state.y[state.npts] = birthprop.v;
	if(marked) 
	  state.marks[state.npts] = birthprop.mrk;
	state.npts     = state.npts + 1;
#ifdef MHDEBUG
	Rprintf("\tnpts=%d\n", state.npts);
#endif
      } else if(itype==DEATH) { 
	/* Death transition */
	/* delete point x[ix], y[ix] */
	if(needupd)
	  (*(cif.update))(state, deathprop, cdata);
	ix = deathprop.ix;
	state.npts = state.npts - 1;
#ifdef MHDEBUG
	Rprintf("implementing death of point %d\n", ix);
	Rprintf("\tnpts=%d\n", state.npts);
#endif
	if(ix < state.npts) {
	  if(!marked) {
	    for(j=ix; j < state.npts; j++) {
	      state.x[j] = state.x[j+1];
	      state.y[j] = state.y[j+1];
	    }
	  } else {
	    for(j = ix; j < state.npts; j++) {
	      state.x[j] = state.x[j+1];
	      state.y[j] = state.y[j+1];
	      state.marks[j] = state.marks[j+1];
	    }
	  }
	}
      } else {              
	/* Shift transition */
	/* Shift (x[ix], y[ix]) to (u,v) */
#ifdef MHDEBUG
	if(marked) 
	  Rprintf("implementing shift from %d = (%lf, %lf)[%d] to (%lf, %lf)[%d]\n", 
		  deathprop.ix, deathprop.u, deathprop.v, deathprop.mrk,
		  shiftprop.u, shiftprop.v, shiftprop.mrk);
	else 
	  Rprintf("implementing shift from %d = (%lf, %lf) to (%lf, %lf)\n", 
		  deathprop.ix, deathprop.u, deathprop.v,
		  shiftprop.u, shiftprop.v);
	Rprintf("\tnpts=%d\n", state.npts);
#endif
	if(needupd)
	  (*(cif.update))(state, shiftprop, cdata);
	ix = shiftprop.ix;
	state.x[ix] = shiftprop.u;
	state.y[ix] = shiftprop.v;
	if(marked) 
	  state.marks[ix] = shiftprop.mrk;
      }
    } 
#ifdef MHDEBUG
    else Rprintf("rejected\n");
#endif
  }

  /* relinquish random number generator */
  PutRNGstate();

  /* return list(x,y) or list(x,y,marks) */
  PROTECT(xout = NEW_NUMERIC(state.npts));
  PROTECT(yout = NEW_NUMERIC(state.npts));
  xx = NUMERIC_POINTER(xout);
  yy = NUMERIC_POINTER(yout);
  for(j = 0; j < state.npts; j++) {
    xx[j] = state.x[j];
    yy[j] = state.y[j];
  }
  if(marked) {
    PROTECT(mout = NEW_INTEGER(state.npts));
    mm = INTEGER_POINTER(mout);
    for(j = 0; j < state.npts; j++) 
      mm[j] = state.marks[j];
  }
  if(!marked) {
    PROTECT(out = NEW_LIST(2));
    SET_VECTOR_ELT(out, 0, xout);
    SET_VECTOR_ELT(out, 1, yout);
    UNPROTECT(19);  /* 16 arguments plus xout, yout, out */
  } else {
    PROTECT(out = NEW_LIST(3)); 
    SET_VECTOR_ELT(out, 0, xout);
    SET_VECTOR_ELT(out, 1, yout); 
    SET_VECTOR_ELT(out, 2, mout);
    UNPROTECT(20);  /* 16 arguments plus xout, yout, mout, out */
  }

  return(out);
}
