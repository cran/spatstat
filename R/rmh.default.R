#
# $Id: rmh.default.R,v 1.18 2004/11/30 19:24:36 adrian Exp adrian $
#
rmh.default <- function(model,start,control=NULL, verbose=TRUE, ...) {
#
# Function rmh.  To simulate realizations of 2-dimensional point
# patterns, given the conditional intensity function of the 
# underlying process, via the Metropolis-Hastings algorithm.
#
# model:   cif par w trend types
# start:   n.start x.start iseed
# control: p q nrep expand periodic ptypes fixall nverb
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

# Name lists --- these will need to be modified when new
# conditional intensity functions are added to the repertoire.
cif.list <- c('strauss','straush','sftcr','straussm','straushm',
               'dgs','diggra','geyer','lookup')
mt.list <- c("straussm","straushm")

# Check that cif is available:
cif <- model$cif
if(!is.loaded(symbol.For(cif)))
	stop(paste("Unrecognized cif: ",cif,".\n",sep=""))

# Turn the name of the cif into a number
nmbr <- match(cif,cif.list)
if(is.na(nmbr)) stop("Name of cif not in cif list.\n")

# Set the 3-vector of integer seeds needed by the subroutine arand:
iseed <- if(is.null(start$iseed)) sample(1:1000000,3) else start$iseed
iseed.save <- iseed

# The R in-house random number/sampling system gets used below.
# Therefore we set an argument for set.seed() --- from iseed --- so
# that if iseed was supplied, then it determines all the randomness,
# including ``proposal points'' and possibly the starting state.
# Thus supplying the same ``iseed'' in a repeat simulation will
# give ***exactly*** the same result.)

	rrr <- .Fortran(
		"arand",
		ix=as.integer(iseed[1]),
		iy=as.integer(iseed[2]),
		iz=as.integer(iseed[3]),
		rand=double(1),
                PACKAGE="spatstat"
	)
	build.seed <- round(rrr$rand*1e6)
	iseed <- unlist(rrr[c("ix","iy","iz")])
	set.seed(build.seed)

# Make sure that the control arguments are adequately specified.
if(is.null(control)) {
	control <- list(
			p=0.9,
			q=0.5,
			nrep=5e5,
			fixall=FALSE,
			periodic=FALSE,
			nverb=0
		   )
}

# Check that precisely one method of specifying the starting
# configuration is provided.  If x.start is given, coerce it into
# a point pattern.  If x.start is specified, model$w may be used to
# specify the simulation window (if this is absent from x.start)
# or a window to which the pattern will be clipped at the finish.
# Dig out these (possibly different) windows.

if(is.null(start$n.start)) {
	if(is.null(start$x.start))
		stop("No starting method provided.\n")
	else {
		whinge <- "Window model$w is not a subwindow of x.start$window.\n"
		use.n <- FALSE
		x.start <- start$x.start
		if(is.null(x.start$window)) {
			if(is.null(model$w))
				stop("No window specified in which to simulate.\n")
			w.sim <- w.clip <- model$w
			x.start <- as.ppp(x.start,w.sim)
		} else {
			x.start <- as.ppp(x.start)
			w.sim   <- x.start$window
			w.clip  <- if(is.null(model$w)) w.sim else model$w
			if(is.null(model$w))
				w.clip <- w.sim
			else {
				if(is.subset.owin(model$w,w.sim)) w.clip <- model$w
				else stop(whinge)
			}
		}
	}
} else {
	if(is.null(start$x.start)) {
		if(is.null(model$w))
			stop("No window specified in which to simulate.\n")
		use.n <- TRUE
		n.start <- start$n.start
		w.clip  <- w.sim <- model$w
	} else stop("Both n.start and x.start were specified.\n")
}

w.sim  <- as.owin(w.sim)
w.clip <- as.owin(w.clip)

# Now (possibly) expand the window; to do this we need to deal
# the expand argument.  Dealing with expand requires that we
# check on several other things.
expand <- control$expand

# Expand must be 1 if there is a trend given by an image.
trendy <- !is.null(model$trend)
if(trendy) {
	trim <- is.im(model$trend) | any(unlist(lapply(model$trend,is.im)))
	if(trim) {
		if(is.null(expand)) expand <- 1
		if(expand > 1)
			stop(paste("When there is trend given as an image,",
                                   "expand must be 1.\n"))
	}
}

# Expand must be 1 if we are conditioning on the number of points.
# cond = 1 <--> no conditioning
# cond = 2 <--> conditioning on n = number of points
# cond = 3 <--> conditioning on the number of points of each type.
p <- control[[match("p",names(control))]]
if(is.null(p)) p <- 0.9
fixall  <- if(is.null(control$fixall)) FALSE else control$fixall
cond    <- 2 - (p<1) + fixall - fixall*(p<1)
bwhinge <- "When conditioning on the number of points,"
if(cond > 1) {
        if(is.null(expand)) expand <- 1
        if(expand > 1)
                stop(paste(bwhinge,"expand must be 1.\n"))
}

# Also the expand argument must be 1 if we are using x.start.
# In a similar vein, if we are using x.start and we are conditioning
# on the number of points, no clipping of the final result
# should be done.  (Nothing to do with expand, as such, but
# we might as well check on this here.)
if(!use.n) {
	if(is.null(expand)) expand <- 1
	else if(expand > 1)
		stop("When using x.start, expand must be 1.\n")
	if(cond > 1 & !identical(all.equal(w.clip,w.sim),TRUE))
		stop(paste(bwhinge,"we cannot clip the result to another window.\n"))
}

# If periodic is TRUE, expand must be 1.  While we're at it,
# if periodic is TRUE check that the window is rectangular and set
# the period.
periodic <- if(is.null(control$periodic)) FALSE else control$periodic
if(periodic) {
	if(is.null(expand)) expand <- 1
        else if(expand > 1)
		stop("Must have expand=1 for periodic simulation.\n")
	if(w.sim$type != "rectangle")
		stop("Need rectangular window for periodic simulation.\n")
	period <- c(w.sim$xrange[2] - w.sim$xrange[1],
                    w.sim$yrange[2] - w.sim$yrange[1])
} else period <- c(-1,-1)

# At this stage we have no reason not to expand the window, and
# if expand has not been specified we let it default to 2.
if(is.null(expand)) expand <- 2

# Now build the expanded window within which to suspend the actual
# window of interest (in order to approximate the simulation of a
# windowed process, rather than a process existing only in the given
# window.  If expand == 1, then we are simulating the latter.  The
# larger ``expand'' is, the better we approximate the former.  Note
# that any value of ``expand'' smaller than 1 is treated as if it
# were 1.
if(expand>1) { # Note that if expand > 1, then use.n is TRUE!
       	rw <- c(w.sim$xrange,w.sim$yrange) # Enclosing box.
	xdim <- rw[2]-rw[1]
	ydim <- rw[4]-rw[3]
	fff  <- (sqrt(expand)-1)*0.5
	ew   <- as.owin(c(rw[1] - fff*xdim,rw[2] + fff*xdim,
                          rw[3] - fff*ydim,rw[4] + fff*ydim))
	n.start <- ceiling(n.start*area.owin(ew)/area.owin(w.sim))
		# The foregoing is incorrect if there is a trend/are trends,
		# but in that case it's just too bloody complicated for it
		# to be worthwhile to try to do the correct thing.
}
else ew <- w.sim

# Determine if the model is multitype.  If so, set up the number
# of types and the types themselves.
if(use.n) {
	mtype <- (!is.null(model$types)) | (length(model$par$beta) > 1)
	npts <- sum(n.start) # Redundant if n.start is scalar; no harm, but.
} else {
	mtype <- is.marked(x.start)
	npts  <- x.start$n
}

if(mtype) {
	if(is.na(match(model$cif, mt.list)))
		stop(paste("Conditional intensity function", cif,
                           "is not multitype.\n"))
	ptypes <- control$ptypes
	if(use.n) {
		types <- model$types
		if (is.null(types)) {
			ntypes <- length(model$par$beta)
			types <- 1:ntypes
		} else {
			ntypes <- length(types)
			if (ntypes != length(model$par$beta))
				stop(paste("Mismatch in lengths of model$types",
                                           "and model$par$beta.\n"))
		}

# If we are conditioning on the number of points of each type, make
# sure that the length of the vector of these numbers is the same as
# ntypes.
		if(cond == 3 & length(n.start) != ntypes)
			stop("Length of n.start not equal to number of types.\n")
# Build ptypes and marks
		if(is.null(ptypes)) ptypes <- rep(1/ntypes,ntypes)
		marks <- if(fixall) rep(1:ntypes,n.start)
			 else sample(1:ntypes,npts,TRUE,ptypes)
	} else {
		marks  <- x.start$marks
                types  <- levels(marks)
                ntypes <- length(types)
                marks  <- match(marks,levels(marks)) # Makes marks into an integer
                                                     # vector for passing to
                                                     # subroutine methas.
		if(is.null(ptypes)) ptypes <- table(marks)/npts
	}

# Check for compatibility of ntypes and length of ptypes.
	if(length(ptypes) != ntypes | sum(ptypes) != 1)
		stop("Argument ptypes is mis-specified.\n")
} else {
	ntypes <- 1
	ptypes <- 1
	marks  <- 0
}

# Warn about a silly value of fixall:
if(fixall & (ntypes==1 | p < 1))
	warning("Setting fixall = TRUE is silly when ntypes = 1 or p < 1.\n")

# Integral of trend over the expanded window (or area of window):
# Iota == Integral Of Trend (or) Area.
if(trendy) {
	if(is.function(model$trend) | is.im(model$trend)) {
		tmp  <- as.im(model$trend,ew)[ew, drop=FALSE]
                sump <- summary(tmp)
		iota <- sump$integral
		tmax <- sump$max
	} else {
		if(is.list(model$trend)) {
			if(length(model$trend) != ntypes)
				stop(paste("Mismatch of ntypes and",
                                           "length of trend list.\n"))
			iota <- numeric(ntypes)
			tmax <- numeric(ntypes)
			for(i in 1:ntypes) {
				tmp <- as.im(model$trend[[i]],ew)[ew, drop=FALSE]
                                sump <- summary(tmp)
				iota[i] <- sump$integral
				tmax[i] <- sump$max
			}
		}
	}
} else {
	iota <- area.owin(ew)
	tmax <- NULL
}

# Build starting state consisting of x, y and possibly marks.
if(use.n) {
	xy   <- if(trendy) rpoint.multi(npts,model$trend,tmax,factor(marks),ew,...)
			else runifpoint(npts,ew,...)
	x <- xy$x
	y <- xy$y
} else {
	x <- x.start$x
	y <- x.start$y
}

# Do some rudimentary pre-processing of the par arguments. Save
# the supplied values of par first.

par.save <- par <- model$par

# 1. Strauss.
if(cif=="strauss") {
	nms <- c("beta","gamma","r")
	if(sum(!is.na(match(names(par),nms))) != 3) {
		cat("For the strauss cif, par must be a named vector\n")
		cat("with components beta, gamma, and r.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters for strauss cif.\n")
	if(par["gamma"] > 1)
		stop("For Strauss processes, gamma must be <= 1.\n")
	par <- par[nms]
}

# 2. Strauss with hardcore.
if(cif=="straush") {
	nms <- c("beta","gamma","r","hc")
	if(sum(!is.na(match(names(par),nms))) != 4) {
		cat("For the straush cif, par must be a named vector\n")
		cat("with components beta, gamma, r, and hc.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters for straush cif.\n")
	par <- par[nms]
}

# 3. Softcore.
if(cif=="sftcr") {
	nms <- c("beta","sigma","kappa")
	if(sum(!is.na(match(names(par),nms))) != 3) {
		cat("For the sftcr cif, par must be a named vector\n")
		cat("with components beta, sigma, and kappa.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative  parameters for sftcr cif.\n")
	if(par["kappa"] > 1)
		stop("For Softcore processes, kappa must be <= 1.\n")
	par <- par[nms]
}

# 4. Marked Strauss.
if(cif=="straussm") {
	nms <- c("beta","gamma","radii")
	if(!is.list(par) || sum(!is.na(match(names(par),nms))) != 3) {
		cat("For the straussm cif, par must be a named list with\n")
		cat("components beta, gamma, and radii.\n")
		stop("Bailing out.\n")
	}
	beta <- par$beta
	if(length(beta) != ntypes)
		stop("Length of beta does not match ntypes.\n")
	gamma <- par$gamma
	if(!is.matrix(gamma) || sum(dim(gamma) == ntypes) != 2)
		stop("Component gamma of par is of wrong shape.\n")
	r <- par$radii
	if(!is.matrix(r) || sum(dim(r) == ntypes) != 2)
		stop("Component r of par is of wrong shape.\n")
	gamma <- t(gamma)[row(gamma)>=col(gamma)]
	r <- t(r)[row(r)>=col(r)]
	par <- c(beta,gamma,r)
	if(any(par[!is.na(par)]<0))
		stop("Negative  parameters for straussm cif.\n")
}

# 5. Marked Strauss with hardcore.
if(cif=="straushm") {
	nms <- c("beta","gamma","iradii","hradii")
	if(!is.list(par) || sum(!is.na(match(names(par),nms))) != 4) {
		cat("For the straushm cif, par must be a named list with\n")
		cat("components beta, gamma, iradii, and hradii.\n")
		stop("Bailing out.\n")
	}
	beta <- par$beta
	if(length(beta) != ntypes)
		stop("Length of beta does not match ntypes.\n")

	gamma <- par$gamma
	if(!is.matrix(gamma) || sum(dim(gamma) == ntypes) != 2)
		stop("Component gamma of par is of wrong shape.\n")
	gamma[is.na(gamma)] <- 1

	iradii <- par$iradii
	if(!is.matrix(iradii) || sum(dim(iradii) == ntypes) != 2)
		stop("Component iradii of par is of wrong shape.\n")
	iradii[is.na(iradii)] <- 0

	hradii <- par$hradii
	if(!is.matrix(hradii) || sum(dim(hradii) == ntypes) != 2)
		stop("Component hradii of par is of wrong shape.\n")
	hradii[is.na(hradii)] <- 0

	gamma <- t(gamma)[row(gamma)>=col(gamma)]
	iradii <- t(iradii)[row(iradii)>=col(iradii)]
	hradii <- t(hradii)[row(hradii)>=col(hradii)]

	par <- c(beta,gamma,iradii,hradii)
	if(any(par[!is.na(par)]<0))
	if(any(is.na(par)) || any(par<0))
		stop("Some parameters negative or illegally set to NA.\n")
}

# 6. Using the interaction function number 1 from Diggle, Gates, and Stibbard).
if(cif=="dgs") {
	nms <- c("beta","rho")
	if(sum(!is.na(match(names(par),nms))) != 2) {
		cat("For the dgs cif, par must be a named vector\n")
		cat("with components beta and rho.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters for dgs cif.\n")
	par <- par[nms]
}

# 7. Using Diggle-Gratton interaction function (function number 2
#    from Diggle, Gates, and Stibbard).
if(cif=="diggra") {
	nms <- c("beta","kappa","delta","rho")
	if(sum(!is.na(match(names(par),nms))) != 4) {
		cat("For the diggra cif, par must be a named vector\n")
		cat("with components beta, kappa, delta, and rho.\n")
		stop("Bailing out.\n")
        }
	if(any(par<0))
		stop("Negative parameters for diggra cif.\n")
	if(par["delta"] >= par["rho"])
		stop("Radius delta must be less than radius rho.\n")
	par <- par[nms]
}

# 8. The Geyer conditional intensity function.
if(cif=="geyer") {
	nms <- c("beta","gamma","r","sat")
	if(sum(!is.na(match(names(par),nms))) != 4) {
		cat("For the geyer cif, par must be a named vector\n")
		cat("with components beta, gamma, r, and sat.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters for geyer cif.\n")
	if(par["sat"] > .Machine$integer.max-100)
		par["sat"] <- .Machine$integer.max-100
	par <- par[nms]
}

# 9. The ``lookup'' device.  This permits simulating, at least
# approximately, ANY pairwise interaction function model
# with isotropic pair interaction (i.e. depending only on distance).
# The pair interaction function is provided as a vector of
# distances and corresponding function values which are used
# as a ``lookup table'' by the Fortran code.

if(cif=="lookup") {
	nms <- c("beta","h")
	if(!is.list(par) || sum(!is.na(match(names(par),nms))) != 2) {
                cat("For the lookup cif, par must be a named list\n")
                cat("with components beta and h (and optionally r).\n")
                stop("Bailing out.\n")
        }
	beta <- par[["beta"]]
	if(beta < 0)
		stop("Negative value of beta for lookup cif.\n")
	h.init <- par[["h"]]
	r <- par[["r"]]
	if(is.null(r)) {
		if(!is.stepfun(h.init))
			stop(paste("For cif=lookup, if component r of",
				   "par is absent then component h must",
				   "be a stepfun object.\n"))
		if(!is.cadlag(h.init))
			stop(paste("The lookup pairwise interaction step",
			     "function must be right continuous,\n",
			     "i.e. built using the default values of",
                             "the \"f\" and \"right\" arguments for stepfun.\n"))
		r     <- knots(h.init)
		h0    <- get("yleft",envir=environment(h.init))
		h     <- h.init(r)
		nlook <- length(r)
		if(!identical(all.equal(h[nlook],1),TRUE))
			stop(paste("The lookup interaction step function",
                                   "must be equal to 1 for \"large\"",
                                   "distances.\n"))
		if(r[1] <= 0)
			stop(paste("The first jump point (knot) of the lookup",
				   "interaction step function must be",
                                   "strictly positive.\n"))
		h <- c(h0,h)
	} else {
		h     <- h.init
		nlook <- length(r)
		if(length(h) != nlook)
			stop("Mismatch of lengths of h and r lookup vectors.\n")
		if(any(is.na(r)))
			stop("Missing values not allowed in r lookup vector.\n")
		if(is.unsorted(r))
			stop("The r lookup vector must be in increasing order.\n")
		if(r[1] <= 0)
			stop(paste("The first entry of the lookup vector r",
                                   "should be strictly positive.\n"))
		h <- c(h,1)
	}
	if(any(h < 0))
		stop(paste("Negative values in the lookup",
                           "pairwise interaction function.\n"))
	if(h[1] > 0 & any(h > 1))
		stop("Lookup pairwise interaction function does not define a valid point process.\n")
	rmax   <- r[nlook]
	r <- c(0,r)
	nlook <- nlook+1
	deltar <- mean(diff(r))
	if(identical(all.equal(diff(r),rep(deltar,nlook-1)),TRUE)) {
		equisp <- 1
		par <- c(beta,nlook,equisp,deltar,rmax,h)
	} else {
		equisp <- 0
		par <- c(beta,nlook,equisp,deltar,rmax,h,r)
	}
}

# If we are simulating a Geyer saturation process we need some
# ``auxilliary information''.
if(nmbr==8) {
	aux <- .Fortran(
		"initaux",
		nmbr=as.integer(nmbr),
		par=as.double(par),
		period=as.double(period),
		x=as.double(x),
		y=as.double(y),
		npts=as.integer(npts),
		aux=integer(npts),
                PACKAGE="spatstat"
	)$aux
	need.aux <- TRUE
} else {
	aux <- 0
	need.aux <- FALSE
}

# The vectors x and y (and perhaps marks) which hold the generated
# process may grow.  We need to allow storage space for them to grow
# in.  Unless we are conditioning on the number of points, we have no
# real idea how big they will grow.  Hence we start off with storage
# space which has twice the length of the ``initial state'', and
# structure things so that the storage space may be incremented
# without losing the ``state'' which has already been generated.

if(cond == 1) {
	x <- c(x,numeric(npts))
	y <- c(y,numeric(npts))
	if(ntypes>1) marks <- c(marks,numeric(npts))
	if(need.aux) aux <- c(aux,numeric(npts))
	nincr <- 2*npts
} else {
	nincr <- npts
}
npmax <- 0
mrep  <- 1

# Set up remaining control parameters:
q     <- if(is.null(control$q)) 0.5 else control$q
nrep  <- if(is.null(control$nrep)) 5e5 else control$nrep
nverb <- if(is.null(control$nverb)) 0 else control$nverb

# If the pattern is multitype, generate the mark proposals.
mprop <- if(ntypes>1)
	sample(1:ntypes,nrep,TRUE,prob=ptypes) else 0
		
# Generate the ``proposal points'' in the expanded window.
xy <- if(trendy) rpoint.multi(nrep,model$trend,tmax,factor(mprop),ew,...) else
			runifpoint(nrep,ew)
xprop <- xy$x
yprop <- xy$y

# The repetition is to allow the storage space to be incremented if
# necessary.
repeat {
	npmax <- npmax + nincr
# Call the Metropolis-Hastings simulator:
	rslt <- .Fortran(
			"methas",
			nmbr=as.integer(nmbr),
			iota=as.double(iota),
			par=as.double(par),
			period=as.double(period),
			xprop=as.double(xprop),
			yprop=as.double(yprop),
			mprop=as.integer(mprop),
			ntypes=as.integer(ntypes),
			ptypes=as.double(ptypes),
			iseed=as.integer(iseed),
			nrep=as.integer(nrep),
			mrep=as.integer(mrep),
			p=as.double(p),
			q=as.double(q),
			npmax=as.integer(npmax),
			nverb=as.integer(nverb),
			x=as.double(x),
			y=as.double(y),
			marks=as.integer(marks),
			aux=as.integer(aux),
			npts=as.integer(npts),
			fixall=as.logical(fixall),
                        PACKAGE="spatstat"
		)

# If npts > npmax we've run out of storage space.  Tack some space
# onto the end of the ``state'' already generated, increase npmax
# correspondingly, and re-call the methas subroutine.  Note that
# mrep is the number of the repetition on which things stopped due
# to lack of storage; so we start again at the ***beginning*** of
# the mrep repetition.
	npts <- rslt$npts
	if(npts <= npmax) break
	if(verbose) {
		cat('Number of points greater than ',npmax,';\n',sep='')
		cat('increasing storage space and continuing.\n')
	}
	x     <- c(rslt$x,numeric(nincr))
	y     <- c(rslt$y,numeric(nincr))
	marks <- if(ntypes>1) c(rslt$marks,numeric(nincr)) else 0
	aux   <- if(need.aux) c(rslt$aux,numeric(nincr)) else 0
	mrep  <- rslt$mrep
	iseed <- rslt$iseed
	npts  <- npts-1
}

x <- rslt$x[1:npts]
y <- rslt$y[1:npts]
if(mtype) {
	marks <- factor(rslt$marks[1:npts],labels=types)
	xxx <- ppp(x=x, y=y, window=as.owin(ew), marks=marks)
} else xxx <- ppp(x=x, y=y, window=as.owin(ew))

# Now clip the pattern to the ``clipping'' window:
xxx <- xxx[,w.clip]

# Append to the result information about how it was generated.
start <- if(use.n) list(n.start=n.start,iseed=iseed.save)
		else list(x.start=x.start,iseed=iseed.save)
attr(xxx, "info") <- list(model=list(cif=cif,par=par.save,trend=model$trend),
                          start=start,
                          control=list(p=p,q=q,nrep=nrep,expand,periodic,
                                       ptypes=ptypes,fixall=fixall))
return(xxx)
}
