#
#       conspire.S
#
#  $Revision: 1.1 $    $Date: 2001/11/05 05:50:40 $
#
#

conspire <- function(indata, formula, subset=NULL, lty=NULL, col=NULL) {

  lhs <- formula[[2]]
  rhs <- formula[[3]]

  # evaluate expression a in data frame b
  evaluate <- function(a,b) {
    if(exists("is.R") && is.R())
      eval(a, envir=b)
    else
      eval(a, local=b)
  }
  
  lhsdata <- evaluate(lhs, indata)
  rhsdata <- evaluate(rhs, indata)
  
  if(is.vector(lhsdata))
    lhsdata <- matrix(lhsdata, ncol=1)

  if(!is.vector(rhsdata))
    stop("rhs of formula seems not to be a vector")

  if(!is.null(subset)) {
    keep <- if(is.character(subset))
		evaluate(parse(text=subset), indata)
            else
                evaluate(subset, indata)
    lhsdata <- lhsdata[keep, , drop=F]
    rhsdata <- rhsdata[keep]
  }

  ylim <- range(lhsdata,na.rm=T)
  xlim <- range(rhsdata,na.rm=T)

  # work out how to label the plot
  yname <- paste(lhs)
  if(length(yname) > 1 && yname[[1]] == "cbind")
    ylab <- paste(yname[-1], collapse=" , ")
  else
    ylab <- as.character(formula)[2]
  
  xlab <- as.character(formula)[3]

  plot(xlim, ylim, type="n", xlab=xlab, ylab=ylab)

  nplots <- ncol(lhsdata)

  if(is.null(lty))
    lty <- 1:nplots
  if(is.null(col))
    col <- rep(1,nplots)
  
  for(i in 1:nplots)
    lines(rhsdata, lhsdata[,i], lty=lty[i], col=col[i])
}
