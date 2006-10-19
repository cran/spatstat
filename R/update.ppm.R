#
#  update.ppm.R
#
#
#  $Revision: 1.16 $    $Date: 2006/10/16 04:47:26 $
#
#
#

update.ppm <- function(object, ..., fixdummy=TRUE, use.internal=NULL) {
  verifyclass(object, "ppm")

  # inspect model
  old <- summary(object, quick="no prediction")$antiquated
  if(old)
    stop("Model was fitted by an old version of spatstat; cannot be updated")
  call <- object$call
  if(!is.call(call))
    stop(paste("Internal error - object$call is not of class",
               sQuote("call")))

  undecided <- is.null(use.internal) || !is.logical(use.internal)
  force.int   <- !undecided && use.internal
  force.ext   <- !undecided && !use.internal
  if(!force.int) {
    # check for validity of format
    badformat <- damaged.ppm(object)
  }
  if(undecided) {
    use.internal <- badformat
    if(badformat)
      message("object format corrupted; repairing it")
  } else if(force.ext && badformat)
    warning("object format corrupted; try update(object, use.internal=TRUE)")
  
  if(use.internal) {
    # reset the main arguments in the call using the internal data
    call$Q <- data.ppm(object)
    namobj <- names(call)
    if("trend" %in% namobj) call$trend <- object$trend
    if("interaction" %in% namobj) call$interaction <- object$interaction
    if("covariates" %in% namobj) call$covariates <- object$covariates
  }
  
  # handle arguments
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
    stop(paste("I am confused: there are two unnamed arguments",
               "of class", sQuote(cname)))
  if(nok == 0) return(0)
  absent <- !(formalname %in% argsgiven)
  if(!absent)
    stop(paste("I am confused: there is an unnamed argument",
               "of class", sQuote(cname), "which conflicts with the",
               "named argument", sQuote(formalname)))
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
    stop(paste("I am confused: there are", nfound,
               "unnamed arguments of different classes (",
               paste(sQuote(cnames(pozzie[found])), collapse=", "),
               ") which could be interpreted as",
               sQuote(formalname)))
}
    

damaged.ppm <- function(object) {
  # guess whether the object format has been damaged
  # e.g. by dump/restore
  gf <- getglmfit(object)
  badfit <- !is.null(gf) && !inherits(gf$terms, "terms")
  if(badfit)
    return(TRUE)
  Qcall <- object$call$Q
  cf <- object$callframe
  if(is.null(cf)) {
    # Old format of ppm objects
    if(is.name(Qcall) && !exists(paste(Qcall)))
      return(TRUE)
    Q <- eval(Qcall)
  } else {
    # New format of ppm objects
    if(is.name(Qcall) && !exists(paste(Qcall), cf))
      return(TRUE)
    Q <- eval(Qcall, cf)
  }
  badQ <- is.null(Q) || !(inherits(Q, "ppp") || inherits(Q,"quad"))
  return(badQ)
}
