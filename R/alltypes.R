#
#      alltypes.R
#
#   $Revision: 1.12 $   $Date: 2006/05/31 04:03:17 $
#
#
alltypes <- function(pp, fun="K", dataname=NULL,verb=FALSE) {
#
# Function 'alltypes' --- calculates a summary function for
# each type, or each pair of types, in a multitype point pattern
#
  verifyclass(pp,"ppp")
  if(!is.character(fun))
    stop("\`fun\' should be a character string")
  
  
# select appropriate statistics
  
  wrong <- function(...) {stop("Internal error!")}
  S <- Si <- Sii <- Sij <- wrong
  
  switch(fun,
         F = {
           indices <- 1
           S  <- Fest
           Si <- function(X, i, ...) { Fest(split(X)[[i]], ...) }
         },
         G = , Gcross = {
           fun <- "Gcross"
           indices <- 2
           S   <- Gest
           Sii <- function(X, i, ...) { Gest(split(X)[[i]], ...) }
           Sij <- function(X, i, j, ...) { Gcross(X, i, j, ...) }
         },
         J = , Jcross = {
           fun <- "Jcross"
           indices <- 2
           S   <- Jest
           Sii <- function(X, i, ...) { Jest(split(X)[[i]], ...) }
           Sij <- function(X, i, j, ...) { Jcross(X, i, j, ...) }
         },
         K = , Kcross = {
           fun <- "Kcross"
           indices <- 2
           S   <- Kest
           Sii <- function(X, i, ...) { Kest(split(X)[[i]], ...) }
           Sij <- function(X, i, j, ...) { Kcross(X, i, j, ...) }
         },
         Gdot = {
           indices <- 1
           S  <- Gest
           Si <- Gdot
         },
         Jdot = {
           indices <- 1
           S  <- Jest
           Si <- Jdot
         },
         Kdot = {
           indices <- 1
           S  <- Kest
           Si <- Kdot
         },
         stop(paste("Unrecognised function name: \`", fun, "\'\n"))
         )

# inspect the possible types  
  if(!is.marked(pp)) {
    um <- 1
    nm <- 1
    indices <- 0
    Tij <- function(X, i, j, ...) { S(X, ...) }
  } else {
    if(!is.factor(pp$marks))
      stop("the marks must be a factor")
    um <- levels(pp$marks)
    nm <- length(um)
    if(indices == 1)
      Tij <- function(X, i, j, ...) { Si(X, i, ...) }
    else
      Tij <- function(X, i, j, ...) {
        if(i == j) Sii(X, i, ...) else Sij(X, i, j, ...)
      }
  }

# build 'fasp' object
  fns  <- list()
  deform <- list()
  
  if(indices <= 1) {
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
            Tij(pp, um[i], um[j], eps=NULL)
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
