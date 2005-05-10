#
#  update.ppm.R
#
#
#  $Revision: 1.8 $    $Date: 2005/04/29 21:22:04 $
#
#
#

update.ppm <- function(object, ..., fixdummy=TRUE) {
  verifyclass(object, "ppm")

  # inspect model
  old <- summary(object, quick="no prediction")$antiquated
  if(old)
    stop("Model was fitted by an old version of spatstat; cannot be updated")
  call <- object$call
  if(!is.call(call))
    stop("Internal error - object$call is not of class \"call\"")

  # arguments
  aargh <- list(...)

  if(length(aargh) == 0)
    return(eval(call, parent.frame()))

  # split named and unnamed arguments
  nama <- names(aargh)
  named <- if(is.null(nama)) rep(FALSE, length(aargh)) else (nama != "")
  namedargs <- aargh[named]
  unnamedargs <- aargh[!named]
  nama <- names(namedargs)
  
  if(any(named)) {
    # any named arguments that were also present in the original call
    # override their original values
    existing <- !is.na(match(nama, names(call)))
    for (a in nama[existing]) call[[a]] <- aargh[[a]]
  
    # add any named arguments not present in the original call
    if (any(!existing)) {
      call <- c(as.list(call), namedargs[!existing])
      call <- as.call(call)
    }
  }
  if(any(!named)) {
    # some objects identified by their class
    if(n <- sp.foundclasses(c("ppp", "quad"), unnamedargs, "Q", nama))
       call$Q <- unnamedargs[[n]]
    if(n<- sp.foundclass("interact", unnamedargs, "interaction", nama))
       call$interaction <- unnamedargs[[n]]
    if(n<- sp.foundclass("formula", unnamedargs, "trend", nama))
       call$trend <- unnamedargs[[n]]
    if(n<- sp.foundclasses(c("data.frame", "im"), unnamedargs, "covariates", nama))
       call$covariates <- unnamedargs[[n]]
  }
  
  # *************************************************************
  # ****** Special action when Q is a point pattern *************
  # *************************************************************
  if(fixdummy && inherits((X <- call$Q), "ppp")) {
    # Instead of allowing default.dummy(X) to occur,
    # explicitly create a quadrature scheme from X,
    # using the same dummy points and weight parameters
    # as were used in the fitted model 
    Qold <- quad.ppm(object)
    Dum <- Qold$dummy
    wpar <- Qold$param$weight
    Qnew <- do.call("quadscheme", append(list(X, Dum), wpar))
    # replace X by new Q
    call$Q <- Qnew
  }

  # finally call ppm
  return(eval(call, parent.frame()))
}

sp.foundclass <- function(cname, inlist, formalname, argsgiven) {
  ok <- unlist(lapply(inlist, inherits, what=cname))
  nok <- sum(ok)
  if(nok > 1)
    stop(paste("I\'m confused: there are two unnamed arguments ",
               "of class \"", cname, "\"", sep=""))
  if(nok == 0) return(0)
  absent <- !(formalname %in% argsgiven)
  if(!absent)
    stop(paste("I\'m confused: there is an unnamed argument ",
               "of class \"", cname, "\" which conflicts with the",
               "named argument \"", formalname, "\"", sep=""))
  theposition <- seq(ok)[ok]
  return(theposition)
}

sp.foundclasses <- function(cnames, inlist, formalname, argsgiven) {
  pozzie <- logical(length(cnames))
  for(i in seq(cnames))
    pozzie[i] <- sp.foundclass(cnames[i],  inlist, formalname, argsgiven)
  found <- (pozzie > 0)
  nfound <- sum(found)
  if(nfound == 0)
    return(0)
  else if(nfound == 1)
    return(pozzie[found])
  else
    stop(paste("I\'m confused: there are ", nfound,
               " unnamed arguments of different classes (\`",
               paste(cnames(pozzie[found]), collapse="\', \`"),
               "\') which could be interpreted as \"",
               formalname, "\"", sep=""))
}
    
