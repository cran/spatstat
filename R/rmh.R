#
#
#    rmh.R
#
#    $Revision: 1.5 $     $Date: 2002/05/13 12:41:10 $
#
#
#

rmh <- function(cif,par,w,ntypes=0,ptypes=NULL,tpar=NULL,n.start,
                expand=NULL,periodic=FALSE,nrep=1e6,p=0.9,q=0.5,
                iseed=NULL,nverb=0) {
#
# Function rmh.  To simulate realizations of 2-dimensional point
# patterns, given the conditional intensity function of the 
# underlying process, via the Metropolis-Hastings algorithm
#

# Check that cif is available:
if(!is.loaded(symbol.For(cif)))
	stop(paste("Unrecognized cif: ",cif,".\n",sep=""))

# Turn the name of the cif into a number
name.list <- c('strauss','straush','sftcr','straussm','straushm',
               'dig1','dig2','geyer')
nmbr <- match(cif,name.list)
if(is.na(nmbr)) stop("Name of cif not in name list.\n")

# Check for compatibility of ntypes and length of ptypes.
if(ntypes <= 1) ptypes <- 1
else {
	if(is.null(ptypes)) ptypes <- rep(1/ntypes,ntypes)
	if(length(ptypes) != ntypes | sum(ptypes) != 1)
		stop("Arguments ntypes and/or ptypes do not make sense.\n")
}

# Turn w into an object of class "owin" if necessary; then
# turn it (back) into a vector of length 4 determining the
# enclosing box (which may simply be the window if the
# original window was rectangular).
w.save <- as.owin(w)
rw <- c(w.save$xrange,w.save$yrange)

# Set the ``period''; if periodic make sure expand is 1.
if(periodic) {
	if(w.save$type != "rectangle") {
		whinge <- paste("\"periodic\" only makes sense",
                                "for rectanglar windows.\n")
		stop(whinge)
	}
	if(is.null(expand)) expand <- 1
	else if(expand > 1) stop("If periodic, expand must be 1.\n")
	period <- c(rw[2] - rw[1], rw[4] - rw[3])
} else {
	if(is.null(expand)) expand <- 2
	period <- c(-1,-1)
}

# Adjust n.start; it is/should be given as the ``expected'' number
# of points in the actual window.  It should be magnified to the
# expected number of points in the bounding box ``rw''.
a1 <- area.owin(w.save)
a2 <- area.owin(as.owin(rw))
n.start <- ceiling(a2*n.start/a1)

# Set need.aux to FALSE (it gets set to TRUE iff cif == 'geyer').
need.aux <- FALSE

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
	need.aux <- TRUE
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

# Set the 3-vector of integer seeds needed by the subroutine arand:
if(is.null(iseed)) iseed <- sample(1:1000000,3)
iseed.save <- iseed

# Prepare vectors x and y (and perhaps marks) to hold the generated
# process; note that we are guessing at how big they will need to be.
# We start off with twice the length of the ``initial state'',
# and structure things so that the storage space may be incremented
# without losing the ``state'' which has already been generated.

npts  <- if(expand > 1) ceiling(expand*n.start) else n.start
n2    <- 2*npts
x     <- numeric(n2)
y     <- numeric(n2)
marks <- if(ntypes > 1) numeric(n2) else 0
aux   <- if(need.aux) numeric(n2) else 0
npmax <- 0
mrep  <- 0 # Indicates starting out; subroutine methas will
           # reset mrep to 1 after generating ``initial state''.

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

# The repetion is to allow the storage space to be incremented if
# necessary.
repeat {
	npmax <- npmax + n2
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
			npts=as.integer(npts)
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
	x     <- c(rslt$x,numeric(n2))
	y     <- c(rslt$y,numeric(n2))
	marks <- if(ntypes>1) c(rslt$marks,numeric(n2)) else 0
	aux   <- if(need.aux) c(rslt$aux,numeric(n2)) else 0
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
                  nrep=nrep,p=p,q=q,expand=expand,periodic=periodic,
                  iseed=iseed.save)
class(tmp) <- "ppp"
tmp
}
