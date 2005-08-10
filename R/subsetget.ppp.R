"[<-.ppp" <-
  function(x, subset, window, value) {
    verifyclass(x, "ppp")
    verifyclass(value, "ppp")
    
    thin <- !missing(subset)
    trim <- !missing(window)
    if(trim)
      verifyclass(window, "owin")

    if(!thin && !trim)
      return(value)

    # determine index subset
    SUB <- seq(x$n)
    
    if(thin) 
      SUB <- SUB[subset]
    if(trim) {
      xsub <- x[SUB]
      xsubok <- inside.owin(xsub$x, xsub$y, window)
      SUB <- SUB[xsubok]
    }

    # anything to replace?
    if(length(SUB) == 0)
      return(x)

    # sanity checks
    if(any(is.na(SUB)))
      stop("Invalid subset: the resulting subscripts include NA\'s")
    
    if(!is.marked.ppp(x) && is.marked.ppp(value))
      warning("The replacement points have marks -- ignored them")
    
    # exact replacement?
    if(value$n == length(SUB)) {
      x$x[SUB] <- value$x
      x$y[SUB] <- value$y
      if(is.marked.ppp(x) && is.marked.ppp(value))
          x$marks[SUB] <- value$marks
    } else {
      if(is.marked.ppp(x) && !is.marked.ppp(value))
        stop("Replacement point pattern must be marked")
      x <- superimpose(x[-SUB], value)
    }
    
    return(x)
}
