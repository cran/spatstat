#
#  update.ppm.R
#
#
#  $Revision: 1.1 $    $Date: 2004/01/27 07:07:02 $
#
#
#

update.ppm <- function(object, ...,
                       Q, trend, interaction, data,
                       correction, rbord, use.gam) {
  verifyclass(object, "ppm")

  aargh <- list(...)

  # Default values for some formal arguments
  defaults <- list(Q=quad.ppm(object),
                   trend=object$trend,
                   interaction=object$interaction,
                   correction=object$correction,
                   rbord=object$rbord)

  matchedargs <- defaults

  # Match named arguments
  # (note: if argument is present and equals NULL, this deletes it from list)
  if(!missing(Q)) matchedargs$Q <- Q
  if(!missing(trend)) matchedargs$trend <- trend
  if(!missing(interaction)) matchedargs$interaction <- interaction
  if(!missing(data)) matchedargs$data <- data
  if(!missing(correction)) matchedargs$correction <- correction
  if(!missing(rbord)) matchedargs$rbord <- rbord
  if(!missing(use.gam)) matchedargs$rbord <- use.gam
  
  # Some formal arguments may be recognised implicitly by their class
  foundclass <- function(cname, inlist, formalname, absent) {
    ok <- unlist(lapply(inlist, inherits, what=cname))
    nok <- sum(ok)
    if(nok > 1)
      stop(paste("I\'m confused: there are two unnamed arguments ",
                 "of class \"", cname, "\"", sep=""))
    if(nok == 0) return(0)
    if(!absent)
      stop(paste("I\'m confused: there is an unnamed argument ",
                 "of class \"", cname, "\" which conflicts with the",
                 "named argument \"", formalname, "\"", sep=""))
    theposition <- seq(ok)[ok]
    return(theposition)
  }
  foundclasses <- function(cnames, inlist, formalname, absent) {
    pozzie <- logical(length(cnames))
    for(i in seq(cnames))
      pozzie[i] <- foundclass(cnames[i],  inlist, formalname, absent)
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

  if(length(aargh) > 0) {
    if(n <- foundclasses(c("ppp", "quad"), aargh, "Q", missing(Q)))
       matchedargs$Q <- aargh[[n]]
    if(n<- foundclass("interact", aargh, "interaction", missing(interaction)))
       matchedargs$interaction <- aargh[[n]]
    if(n<- foundclass("formula", aargh, "trend", missing(trend)))
       matchedargs$trend <- aargh[[n]]
    if(n<- foundclass("data.frame", aargh, "data", missing(data)))
       matchedargs$data <- aargh[[n]]
  }
  
  # *************************************************************
  # ****** Special action when Q is a point pattern *************
  # *************************************************************
  if(!is.null(X <- matchedargs$Q) && inherits(X, "ppp")) {
    # Instead of allowing default.dummy(X) to occur,
    # explicitly create a quadrature scheme from X,
    # using the same dummy points and weight parameters
    # as were used in the fitted model 
    Qold <- quad.ppm(object)
    Dum <- Qold$dummy
    wpar <- Qold$param$weight
    Qnew <- do.call("quadscheme", append(list(X, Dum), wpar))
    # replace X by new Q
    matchedargs$Q <- Qnew
  }

  # finally call mpl
  result <- do.call("mpl", matchedargs)
  return(result)
}
