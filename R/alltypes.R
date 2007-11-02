#
#      alltypes.R
#
#   $Revision: 1.15 $   $Date: 2007/10/30 17:15:35 $
#
#
alltypes <- function(pp, fun="K", dataname=NULL,verb=FALSE) {
#
# Function 'alltypes' --- calculates a summary function for
# each type, or each pair of types, in a multitype point pattern
#
  verifyclass(pp,"ppp")
  if(!is.character(fun))
    stop(paste(sQuote("fun"), "should be a character string"))
  
  if(is.null(dataname)) dataname <- deparse(substitute(pp))
  
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
         stop(paste("Unrecognised function name:", sQuote(fun)))
         )

# inspect the possible types  
  if(!is.marked(pp)) {
    um <- 1
    nm <- 1
    indices <- 0
    Tij <- function(X, i, j, ...) { S(X, ...) }
  } else {
    ma <- marks(pp)
    if(!is.factor(ma))
      stop("the marks must be a factor")
    um <- levels(ma)
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

  marklabels <- paste(um)
  witch <-
    if(indices == 0) 
      matrix(1, nrow=1, ncol=1, dimnames=list("", ""))
    else if(indices == 1) 
      matrix(1:nm, nrow=nm, ncol=1,
             dimnames=list(marklabels, ""))
    else 
      matrix(1:(nm^2),ncol=nm,nrow=nm, byrow=TRUE,
             dimnames <- list(marklabels, marklabels))

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
  if(nm > 1)
	title <- paste("Array of ",fun," functions for ",
              	dataname,".",sep="")
  else
	title <- paste(fun," function for ",dataname,".",sep="")

  rslt <- fasp(fns, which=witch,
               formulae=deform,
               dataname=dataname,
               title=title)
  return(rslt)
}
