#
#      alltypes.R
#
#   $Revision: 1.9 $   $Date: 2004/01/13 10:05:46 $
#
#
alltypes <- function(pp, fun="K", dataname=NULL,verb=FALSE) {
#
# Function 'alltypes' --- calculates a summary function for
# each type, or each pair of types, in a multitype point pattern
#
  verifyclass(pp,"ppp")

# validate 'fun'
  switch(fun, F={}, G={}, J={}, K={},
         stop("Unrecognized function name: ",fun,".\n"))
  wrong <- function(...) {stop("Internal error!")}

# list all possible types  
  if(!is.marked(pp)) {
    um <- 1
    nm <- 1
  } else {
    if(!is.factor(pp$marks))
      stop("the marks must be a factor")
    um <- levels(pp$marks)
    nm <- length(um)
  }
  
# get sensible 'r' values
  brks <- handle.r.b.args(r=NULL, breaks=NULL, pp$window, eps=NULL)
  r    <- brks$r

  # select appropriate statistics
  F1 <- switch(fun,F=Fest,G=Gest,J=Jest,K=Kest, wrong)
  F2 <- switch(fun,F={},G=Gcross,J=Jcross,K=Kcross, wrong)

# build 'fasp' object
  fns  <- list()
  deform <- list()
  
  if(fun=="F") {
    witch <- matrix(1:nm,ncol=1,nrow=nm)
    names(witch) <- um
    titles <- if(nm > 1) as.list(paste("mark =", um)) else list("")
  } else {
    witch <- matrix(1:(nm^2),ncol=nm,nrow=nm,byrow=TRUE)
    dimnames(witch) <- list(um, um)
    titles <- if(nm > 1)
      as.list(paste("(", um[t(row(witch))], ",", um[t(col(witch))], ")", sep=""))
    else
      list("")
  }

  # compute function array
  k   <- 0

  for(i in 1:nrow(witch)) {
	for(j in 1:ncol(witch)) {
          if(verb) cat("i =",i,"j =",j,"\n")
          k <- k+1
          fns[[k]] <- currentfv <- 
            if(nm == 1) # univariate pattern
              F1(pp,eps=NULL,breaks=brks)
            else if(fun=="F" | i==j) # F_i or G_ii, J_ii, K_ii
              F1(pp[pp$marks==um[i]], eps=NULL,breaks=brks)
            else 
              F2(pp,um[i],um[j], eps=NULL, breaks=brks)
          deform[[k]] <- attr(currentfv, "fmla")
        }
      }

  # wrap up into 'fasp' object
  if(is.null(dataname)) dataname <- deparse(substitute(pp))

  if(nm > 1)
	title <- paste("Array of ",fun," functions for ",
              	dataname,".",sep="")
  else
	title <- paste(fun," function for ",dataname,".",sep="")

  rslt <- fasp(fns, titles, deform, witch, dataname, title)
  return(rslt)
}
