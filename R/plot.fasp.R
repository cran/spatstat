#
#   plot.fasp.R
#
#   $Revision: 1.5 $   $Date: 2002/01/18 06:27:22 $
#
plot.fasp <- function(x,formula=NULL,subset=NULL,lty=NULL,
                      col=NULL,title=NULL,...) {

# If the formula is null, look for a default formula in x:
	if(is.null(formula)) {
		if(is.null(x$default.formula))
			stop("No formula supplied.\n")
		formula <- x$default.formula
	}
# The formula should be a single formula or a list of formulae.
# If it is a single formula, wrap it up in a list so that all
# formula arguments can be treated consistently.
        if(!is.list(formula)) formula <- list(formula)

# Check on the length of the formula argument.
nf <- length(formula)
if(nf > 1) {
	if(nf != length(x$fns))
		stop("Wrong number of entries in formula argument.\n")
	mfor <- T
} else mfor <- F

# Check on the length of the subset argument.
ns <- length(subset)
if(ns > 1) {
	if(ns != length(x$fns))
		stop("Wrong number of entries in subset argument.\n")
	msub <- T
} else msub <- F

# Set up the array of plotting regions.
	mfrow.save <- par("mfrow")
	oma.save   <- par("oma")
	on.exit(par(mfrow=mfrow.save,oma=oma.save))
	which <- x$which
	m  <- nrow(which)
	n  <- ncol(which)
	nm <- n * m
	par(mfrow=c(m,n))
        # decide whether panels require subtitles
        subtit <- (nm > 1) || !(is.null(x$titles[[1]]) || x$titles[[1]] == "")
	if(nm>1) par(oma=c(0,3,4,0))
        else if(subtit) par(oma=c(3,3,4,0))
        
# Run through the components of the structure x, plotting each
# in the appropriate region, according to the formula.
	k <- 0
	for(i in 1:m) {
		for(j in 1:n) {
# Now do the actual plotting.
			k <- which[i,j]
			if(is.na(k)) plot(0,0,type='n',xlim=c(0,1),
					  ylim=c(0,1),axes=F,xlab='',ylab='')
			else {
				dat.loc <- x$fns[[k]]
				sij <- if(msub) subset[[k]] else subset
				fij <- if(mfor) formula[[k]] else formula[[1]]
				conspire(dat.loc,fij,sij,lty,col)

# Add the (sub)title of each plot.
				if(!is.null(x$titles[[k]]))
					title(main=x$titles[[k]])
			}
		}
	}

# Add an overall title.
	if(!is.null(title)) overall <- title
	else if(!is.null(x$title)) overall <- x$title
	else {
		if(nm > 1)
			overall <- "Array of diagnostic functions"
		else
			overall <- "Diagnostic function"
		if(is.null(x$dataname)) overall <- paste(overall,".",sep="")
		else overall <- paste(overall," for ",x$dataname,".",sep="")
		
	}
	if(nm > 1 || subtit)
          mtext(side=3,outer=T,line=1,text=overall,cex=1.2)
	else title(main=overall)
	invisible()
}
