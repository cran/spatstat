#
#   plot.fasp.R
#
#   $Revision: 1.12 $   $Date: 2005/05/06 03:32:43 $
#
plot.fasp <- function(x, formule=NULL, subset=NULL,
                      lty=NULL, col=NULL, title=NULL, ..., samex=TRUE) {

# If no formula is given, look for a default formula in x:
  defaultplot <- is.null(formule)
  if(defaultplot) {
    if(is.null(x$default.formula))
      stop("No formula supplied.\n")
    formule <- x$default.formula
  }

# ensure formulae are given as character strings.  
  formule <- FormatFaspFormulae(formule, "formule")

# Number of formulae should match number of functions.
  nf <- length(formule)
  nfun <- length(x$fns)
  if(nf == 1 && nfun > 1)
    formule <- rep(formule, nfun)
  else if(nf != nfun)
    stop("Wrong number of entries in \`formule\'.\n")

# Check on the length of the subset argument.
  ns <- length(subset)
  if(ns > 1) {
    if(ns != length(x$fns))
      stop("Wrong number of entries in subset argument.\n")
    msub <- TRUE
  } else msub <- FALSE

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

# compute common x axis limits for all plots (in default case)
  if(defaultplot && samex) {
    ends <- lapply(x$fns, function(z) { attr(z, "alim")})
    isnul <- unlist(lapply(ends, is.null))
    ends <- ends[!isnul]
    lo <- max(unlist(lapply(ends, min)))
    hi <- min(unlist(lapply(ends, max)))
    xlim <- c(lo,hi)
  } else xlim <- NULL
          
# Run through the components of the structure x, plotting each
# in the appropriate region, according to the formula.
  k <- 0
  for(i in 1:m) {
    for(j in 1:n) {
# Now do the actual plotting.
      k <- which[i,j]
      if(is.na(k)) plot(0,0,type='n',xlim=c(0,1),
                        ylim=c(0,1),axes=FALSE,xlab='',ylab='', ...)
      else {
        fun <- as.fv(x$fns[[k]])
        fmla <- formule[k] 
        sub <- if(msub) subset[[k]] else subset
        if(is.null(xlim))
          plot(fun, fmla, sub, lty,col, ...)
        else
          plot(fun, fmla, sub, lty,col, xlim=xlim, ...)

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
    mtext(side=3,outer=TRUE,line=1,text=overall,cex=1.2)
  else title(main=overall)

  invisible()
}
