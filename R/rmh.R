rmh <- function(cif,par,w=NULL,ntypes=1,ptypes=NULL,tpar=NULL,n.start=NULL,
                x.start=NULL,fixall=FALSE,expand=NULL,periodic=FALSE,
                nrep=1e6,p=0.9,q=0.5,iseed=NULL,nverb=0) {
#
# Function rmh.  To simulate realizations of 2-dimensional point
# patterns, given the conditional intensity function of the 
# underlying process, via the Metropolis-Hastings algorithm.
#

# Check that cif is available:
if(!is.loaded(symbol.For(cif)))
	stop(paste("Unrecognized cif: ",cif,".\n",sep=""))

# Turn the name of the cif into a number
name.list <- c('strauss','straush','sftcr','straussm','straushm',
               'dig1','dig2','geyer')
nmbr <- match(cif,name.list)
if(is.na(nmbr)) stop("Name of cif not in name list.\n")

# Check that precisely one method of specifying the starting
# configuration is provided.
start.count <- sum(c(is.null(n.start),is.null(x.start)))
if(start.count != 1)
	stop("Specify precisely one of n.start or x.start.")

# Save the name of x.start (may be NULL).
x.start.name <- deparse(substitute(x.start))

# If x.start is given, coerce it into a point pattern.
if(!is.null(x.start)) {
	if(is.null(x.start$w) & !is.null(w))
		x.start <- as.ppp(x.start,w)
	else x.start <- as.ppp(x.start)
}

# Check that a window in which to simulate is provided.
w.count <- sum(c(is.null(w),is.null(x.start)))
if(w.count == 2)
	stop("No window specified in which to simulate.\n")

# Save the window to which the simulated process will be
# clipped at the finish.
w.save <- if(is.null(w)) x.start$w else as.owin(w)

# Determine if the window is rectangular
type    <- if(is.null(x.start)) as.owin(w)$type else x.start$w$type
rectwin <- type == "rectangle"

# What sort of conditioning are we doing?
# cond = 1 <--> no conditioning
# cond = 2 <--> conditioning on n = number of points
# cond = 3 <--> conditioning on the number of points of each type.
cond <- 2 - (p<1) + fixall - fixall*(p<1)

# Now check that arguments make sense for the sort of conditioning
# that we are doing.
# If cond > 1, then
#	o expand must equal 1 (this includes w == x.start$w)
#	o the window must be rectangular
# If cond = 3, then
#	o length(n.start) must equal ntypes
bwinge <- "When conditioning on the number of points,"
if(cond > 1) {
	if(!is.null(expand) && expand > 1)
		stop(paste(bwinge,"expand must be 1.\n"))
	if(w.count == 0 && !identical(all.equal(as.owin(w),x.start$w),TRUE))
		stop(paste(bwinge,"w must be the same as x.start$w.\n"))
	if(!rectwin)
		stop(paste(bwinge,"window must be rectangular.\n"))
}

# Set the 3-vector of integer seeds needed by the subroutine arand:
if(is.null(iseed)) iseed <- sample(1:1000000,3)
iseed.save <- iseed

# If n.start is specified, check that other arguments make
# sense, and construct the starting configuration.
if(!is.null(n.start)) {

# Turn w into a vector of length 4 determining the
# enclosing box (which may simply be the window if the
# original window was rectangular).
        rw <- c(w.save$xrange,w.save$yrange)

# Set the value of expand.
	if(is.null(expand)) expand <- if(cond > 1 | periodic) 1 else 2

# Build the expanded window within which to suspend the actual
# window of interest (in order to approximate the simulation of a
# windowed process, rather than a process existing only in the given
# window.  If expand == 1, then we are simulating the latter.  The
# larger ``expand'' is, the better we approximate the former.  Note
# that any value of ``expand'' smaller than 1 is treated as if it
# were 1.
	if(expand>1) {
		xdim <- rw[2]-rw[1]
		ydim <- rw[4]-rw[3]
		fff  <- (sqrt(expand)-1)*0.5
		erw  <- c(rw[1] - fff*xdim,rw[2] + fff*xdim,
                          rw[3] - fff*ydim,rw[4] + fff*ydim)
	}
	else erw <- rw

# If we are conditioning on the number of points of each type, make
# sure that the length of the vector of these numbers is the same as
# ntypes.
if(cond == 3 & length(n.start) != ntypes)
	stop("Length of n.start not equal to ntypes.\n")

# Even if we are conditioning on the number of points of each type,
# rather than just on the total, we still need the total.  Get the
# total number of points in the starting configuration (redundant
# unless we are conditioning on the number of points of each type).
	npts <- sum(n.start)

# Adjust npts; it is/should be given as the ``expected'' number
# of points in the actual window.  It should be magnified to the
# expected number of points in the bounding box ``rw''.
	a1 <- area.owin(w.save)
	a2 <- area.owin(as.owin(rw))
	npts <- ceiling(a2*npts/a1)

# Check for compatibility of ntypes and length of ptypes.
	if(ntypes <= 1) ptypes <- 1
	else {
		if(is.null(ptypes)) ptypes <- rep(1/ntypes,ntypes)
		if(length(ptypes) != ntypes | sum(ptypes) != 1)
			stop("Problem with arguments ntypes/ptypes.\n")
	}

# Set the random number generation seed, in a manner coordinated
# with the seeds which will be passed on to the subroutine arand:
rrr <- .Fortran(
		"arand",
		ix=as.integer(iseed[1]),
		iy=as.integer(iseed[2]),
		iz=as.integer(iseed[3]),
		rand=double(1)
	)
build.seed <- round(rrr$rand*1e6)
iseed <- unlist(rrr[c("ix","iy","iz")])
set.seed(build.seed)

# Build starting state consisting of x, y and possibly marks.
	npts  <- if(expand > 1) ceiling(expand*npts) else npts
	x <- runif(npts,erw[1],erw[2])
	y <- runif(npts,erw[3],erw[4])
	if(ntypes > 1) {
		marks <- if(fixall) rep(1:ntypes,n.start) else
				    sample(1:ntypes,npts,TRUE,ptypes)
	} else marks <- 0
}

# If x.start is specified, set up other arguments from x.start:
else {
	erw <- x.start$w
	erw <- c(erw$xrange,erw$yrange)
	npts <- x.start$n
	x    <- x.start$x
	y    <- x.start$y
	expand <- 1
	if(is.marked(x.start)) {
		marks  <- x.start$marks
		ntypes <- length(levels(marks))
		marks  <- match(marks,levels(marks))
		if(is.null(ptypes)) ptypes <- table(marks)/npts
	} else {
		ntypes <- 1
		marks  <- 0
		ptypes <- 1
	}
}

# Warn about a silly value of fixall:
if(fixall & (ntypes==1 | p < 1))
	warning("Setting fixall = TRUE is silly when ntypes = 1 or p < 1.\n")

# Set the ``period''.
if(periodic) {
	if(!rectwin)
		stop("Need rectangular window for periodic simulation.\n")
	if(expand>1)
		stop("Must have expand=1 for periodic simulation.\n")
	period <- c(erw[2] - erw[1], erw[4] - erw[3])
} else period <- c(-1,-1)

# Do some rudimentary pre-processing of the par and tpar arguments.
# Save the supplied values of par and tpar first.
par.save  <- par
tpar.save <- tpar

# 1. Strauss.
if(cif=="strauss") {
	if(length(par) != 3) {
		cat("For strauss cif, par should be a vector of\n")
		cat("length 3, consisting of beta, gamma and r.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters.\n")
	if(par[2] > 1)
		stop("For Strauss processes, gamma must be <= 1.\n")
}

# 2. Strauss with hardcore.
if(cif=="straush") {
	if(length(par) != 4) {
		cat("For straush cif, par should be a vector of\n")
		cat("length 4, consisting of beta, gamma, r, and r_0.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters.\n")
}

# 3. Softcore.
if(cif=="sftcr") {
	if(length(par) != 3) {
		cat("For sftcr cif, par should be a vector of\n")
		cat("length 3, consisting of beta, sigma and kappa.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative  parameters.\n")
	if(par[3] > 1)
		stop("For Softcore processes, kappa must be <= 1.\n")
}

# 4. Marked Strauss.
if(cif=="straussm") {
	if(ntypes<=1)
		stop("Argument ntypes must be at least 2 for straussm.\n")
	if(!is.atomic(par)) {
		if(!is.list(par)) {
			cat("For straussm cif, par must be either a\n")
			cat("vector, or a list with components beta,\n")
			cat("gamma, and r.\n")
			stop("Bailing out.\n")
		}
		beta <- par$beta
		if(length(beta) != ntypes)
			stop("Component beta of par is of wrong length.\n")
		gamma <- par$gamma
		if(!is.matrix(gamma) | length(gamma) != ntypes^2)
			stop("Component gamma of par is of wrong shape.\n")
		r <- par$r
		if(!is.matrix(r) | length(r) != ntypes^2)
			stop("Component r of par is the wrong shape.\n")
		gamma <- t(gamma)[row(gamma)>=col(gamma)]
		r <- t(r)[row(r)>=col(r)]
		par <- c(beta,gamma,r)
	}
	if(any(par[!is.na(par)]<0))
		stop("Negative  parameters.\n")
}

# 5. Marked Strauss with hardcore.
if(cif=="straushm") {
	if(ntypes<=1)
		stop("Argument ntypes must be at least 2 for straushm.\n")
	if(!is.atomic(par)) {
		if(!is.list(par)) {
			cat("For straushm cif, par must be either a\n")
			cat("vector, or a list with components beta,\n")
			cat("gamma, r, and rhc.\n")
			stop("Bailing out.\n")
		}
		beta <- par$beta
		if(length(beta) != ntypes)
			stop("Component beta of par is of wrong length.\n")

		gamma <- par$gamma
		if(!is.matrix(gamma) | length(gamma) != ntypes^2)
			stop("Component gamma of par is of wrong shape.\n")

		r <- par$r
		if(!is.matrix(r) | length(r) != ntypes^2)
			stop("Component r of par is the wrong shape.\n")

		rhc <- par$rhc
		if(!is.matrix(rhc) | length(rhc) != ntypes^2)
			stop("Component rhc of par is the wrong shape.\n")

		gamma <- t(gamma)[row(gamma)>=col(gamma)]
		r <- t(r)[row(r)>=col(r)]
		rhc <- t(rhc)[row(rhc)>=col(rhc)]

		par <- c(beta,gamma,r,rhc)
	}
	if(any(par[!is.na(par)]<0))
		stop("Negative  parameters.\n")
}

# 6. Using interaction function number 1 from Diggle, Gates,
#    and Stibbard.

if(cif=="dig1") {
	if(length(par) != 2) {
		cat("For dig1 cif, par should be a vector of\n")
		cat("length 2, consisting of beta and rho.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters.\n")
}

# 7. Using interaction function number 2 from Diggle, Gates,
#    and Stibbard.

if(cif=="dig2") {
	if(length(par) != 4) {
		cat("For dig2 cif, par should be a vector of\n")
		cat("length 4, consisting of beta, kappa, delta\n")
		cat("and rho.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters.\n")
	if(par[3] >= par[4])
		stop("Radius delta must be less than radius rho.\n")
}

# 8. The Geyer conditional intensity function.

if(cif=="geyer") {
	if(length(par) != 4) {
		cat("For the geyer cif, par should be a vector of\n")
		cat("length 4, consisting of beta, gamma, r\n")
		cat("and s.\n")
		stop("Bailing out.\n")
	}
	if(any(par<0))
		stop("Negative parameters.\n")
	if(par[4] > .Machine$integer.max-100)
		par[4] <- .Machine$integer.max-100
}

# Calculate the degree(s) of the polynomial(s) in the log polynomial
# trend(s), if any, and tack it/them onto tpar.  Note that
# if tpar is NULL, this will make tpar into a scalar equal to 0.
if(ntypes > 1 & (!is.null(tpar))) {
	if(!is.atomic(tpar)) {
		if(!is.list(tpar) || length(tpar) != ntypes)
			stop("Argument tpar has wrong form.\n")
		tmp <- NULL
		for(cc in tpar) {
			nc <- length(cc)
			nd <- (sqrt(9+8*nc) - 3)/2
			if(nd%%1 > .Machine$double.eps) {
				stpms <- "Length of a component"
				stpms <- paste(stpms,"does not make sense.\n")
				stop(stpms)
			}
			tmp <- c(tmp,nd)
		}
		tpar <- c(ntypes,tmp,unlist(tpar))
	}
}
else {
	nt <- length(tpar)
	nd <- (sqrt(9+8*nt) - 3)/2
	if(nd%%1 > .Machine$double.eps)
	stop("Length of tpar does not make sense.\n")
	tpar <- c(nd,tpar)
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
		aux=integer(npts)
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

# The repetion is to allow the storage space to be incremented if
# necessary.
repeat {
	npmax <- npmax + nincr
# Call the Metropolis-Hastings simulator:
	rslt <- .Fortran(
			"methas",
			nmbr=as.integer(nmbr),
			rw=as.double(erw),
			par=as.double(par),
			period=as.double(period),
			tpar=as.double(tpar),
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
			fixall=as.logical(fixall)
		)

# If npts > npmax we've run out of storage space.  Tack some space
# onto the end of the ``state'' already generated, increase npmax
# correspondingly, and re-call the hasmet subroutine.  Note that
# mrep is the number of the repetion on which things stopped due
# to lack of storage; so we start again at the ***beginning*** of
# the mrep repetion.
	npts <- rslt$npts
	if(npts <= npmax) break
	cat('Number of points greater than ',npmax,';\n',sep='')
	cat('increasing storage space and continuing.\n')
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
if(ntypes>1) marks <- rslt$marks[1:npts]
tmp <- as.ppp(list(x=x,y=y),W=as.owin(erw))
if(ntypes>1) tmp$marks <- factor(marks)

# Now window the returned process by the original window:
tmp <- tmp[,w.save]

# Append to the result information about how it was generated.
tmp$info <- list(cif=cif,par=par.save,tpar=tpar.save,n.start=n.start,
		  x.start.name=x.start.name,nrep=nrep,p=p,q=q,
                  expand=expand,periodic=periodic,iseed=iseed.save)

class(tmp) <- "ppp"
tmp
}
