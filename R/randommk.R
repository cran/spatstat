#
#
#   randommk.R
#
#   Random generators for MULTITYPE point processes
#
#   $Revision: 1.8 $   $Date: 2005/03/08 21:29:51 $
#
#   rmpoispp()   random marked Poisson pp
#   rmpoint()    n independent random marked points
#   rmpoint.I.allim()  ... internal
#   rpoint.multi()   temporary wrapper 
#
"rmpoispp" <-
  function(lambda, lmax=NULL, win = owin(c(0,1),c(0,1)),
           types, ...) {
    # arguments:
    #     lambda  intensity:
    #                constant, function(x,y,m,...), image,
    #                vector, list of function(x,y,...) or list of images
    #
    #     lmax     maximum possible value of lambda
    #                constant, vector, or list
    #
    #     win     default observation window (of class 'owin')
    #
    #     types    possible types for multitype pattern
    #    
    #     ...     extra arguments passed to lambda()
    #

    # Validate arguments
    is.numvector <- function(x) {is.numeric(x) && is.vector(x)}
    is.constant <- function(x) {is.numvector(x) && length(x) == 1}
    checkone <- function(x) {
      if(is.constant(x)) {
        if(x >= 0) return(TRUE) else stop("Intensity is negative!")
      }
      return(is.function(x) || is.im(x))
    }
    single.arg <- checkone(lambda)
    vector.arg <- !single.arg && is.numvector(lambda) 
    list.arg <- !single.arg && is.list(lambda)
    if(! (single.arg || vector.arg || list.arg))
      stop("argument \'lambda\' not understood")
    
    if(list.arg && !all(unlist(lapply(lambda, checkone))))
      stop("Each entry in the list \'lambda\' must be either a constant, a function or an image")
    if(vector.arg && any(lambda < 0))
      stop("Some entries in the vector \'lambda\' are negative")


    # Determine & validate the set of possible types
    if(missing(types)) {
      if(single.arg)
        stop("\'types\' must be given explicitly\
 when \'lambda\' is a constant, a function or an image")
      else
        types <- seq(lambda)
    } 

    ntypes <- length(types)
    if(!single.arg && (length(lambda) != ntypes))
      stop("The lengths of \'lambda\' and \'types\' do not match")

    factortype <- factor(types, levels=types)

    # Validate `lmax'
    if(! (is.null(lmax) || is.numvector(lmax) || is.list(lmax) ))
      stop("\'lmax\' should be a constant, a vector, a list or NULL")
       
    # coerce lmax to a vector, to save confusion
    if(is.null(lmax))
      maxes <- rep(NULL, ntypes)
    else if(is.numvector(lmax) && length(lmax) == 1)
      maxes <- rep(lmax, ntypes)
    else if(length(lmax) != ntypes)
      stop("The length of \'lmax\' does not match the number of possible types")
    else if(is.list(lmax))
      maxes <- unlist(lmax)
    else maxes <- lmax

    # coerce lambda to a list, to save confusion
    lam <- if(single.arg) lapply(1:ntypes, function(x, y){y}, y=lambda)
           else if(vector.arg) as.list(lambda) else lambda
    
    # Simulate
    for(i in 1:ntypes) {
      if(single.arg && is.function(lambda))
        # call f(x,y,m, ...)
        Y <- rpoispp(lambda, lmax=maxes[i], win=win, types[i], ...)
      else
        # call f(x,y, ...) or use other formats
        Y <- rpoispp(lam[[i]], lmax=maxes[i], win=win, ...)
      Y <- Y %mark% factortype[i]
      X <- if(i == 1) Y else superimpose(X, Y)
    }

    # Randomly permute, just in case the order is important
    permu <- sample(X$n)
    return(X[permu])
}

# ------------------------------------------------------------------------

"rmpoint" <- function(n, f=1, fmax=NULL, 
                      win = unit.square(), 
                      types, ptypes, ...,
                      giveup = 1000, verbose = FALSE) {
  if(!is.numeric(n))
    stop("n must be a scalar or vector")

  if(sum(n) == 0)
    return(runifpoint(0, f=f, fmax=fmax, win=win) %mark% sample(types, 0))

  #############
  
  Model <- if(length(n) == 1) {
    if(missing(ptypes)) "I" else "II"
  } else "III"
  
  ##############  Validate f argument
  is.numvector <- function(x) {is.numeric(x) && is.vector(x)}
  is.constant <- function(x) {is.numvector(x) && length(x) == 1}
  checkone <- function(x) {
    if(is.constant(x)) {
      if(x >= 0) return(TRUE) else stop("Intensity is negative!")
    }
    return(is.function(x) || is.im(x))
  }

  single.arg <- checkone(f)
  vector.arg <- !single.arg && is.numvector(f) 
  list.arg <- !single.arg && is.list(f)
  if(! (single.arg || vector.arg || list.arg))
    stop("argument \'f\' not understood")
    
  if(list.arg && !all(unlist(lapply(f, checkone))))
    stop("Each entry in the list \'f\' must be either a constant, a function or an image")
  if(vector.arg && any(f < 0))
    stop("Some entries in the vector \'f\' are negative")


  ################   Determine & validate the set of possible types
  if(missing(types)) {
    if(single.arg && length(n) == 1)
      stop("\'types\' must be given explicitly\
 when \'f\' is a single number, a function or an image\
 and \`n\' is a single number")
    else if(single.arg)
      types <- seq(n)
    else 
      types <- seq(f)
  }

  ntypes <- length(types)
  if(!single.arg && (length(f) != ntypes))
    stop("The lengths of \'f\' and \'types\' do not match")
  if(length(n) > 1 && ntypes != length(n))
    stop("The lengths of \'n\' and \'types\' do not match")

  factortype <- factor(types, levels=types)
  
  #######################  Validate `fmax'
  if(! (is.null(fmax) || is.numvector(fmax) || is.list(fmax) ))
    stop("\'fmax\' should be a constant, a vector, a list or NULL")
       
  # coerce fmax to a vector, to save confusion
  if(is.null(fmax))
    maxes <- rep(NULL, ntypes)
  else if(is.constant(fmax))
    maxes <- rep(fmax, ntypes)
  else if(length(fmax) != ntypes)
    stop("The length of \'fmax\' does not match the number of possible types")
  else if(is.list(fmax))
    maxes <- unlist(fmax)
  else maxes <- fmax

  # coerce f to a list, to save confusion
  flist <- if(single.arg) lapply(1:ntypes, function(i, f){f}, f=f)
         else if(vector.arg) as.list(f) else f

  #################### START ##################################

  ## special algorithm for Model I when all f[[i]] are images

  if(Model == "I" && all(unlist(lapply(f, is.im))))
    return(rmpoint.I.allim(n, f, types))

  ## otherwise, first select types, then locations given types
  
  if(Model == "I") {
    # Compute approximate marginal distribution of type
    integratexy <- function(f, win, ...) {
      imag <- as.im(f, win, ...)
      summ <- summary(imag)
      summ$integral
    }
    integratexyi <- function(i, f, win, ...) { integratexy(f, win, i, ...) }
    
    fintegrals <- if(vector.arg) f * area.owin(win) else
       if(list.arg) unlist(lapply(flist, integratexy, win=win, ...)) else
       unlist(lapply(1:ntypes, integratexyi, f=f, win=win, ...))
    ptypes <- fintegrals/sum(fintegrals)
  }

  # Generate number of points of each type

  if(Model == "I" || Model == "II") {
    randomtypes <- sample(types, n, prob=ptypes, replace=TRUE)
    nn <- table(randomtypes)
  } else
    nn <- n

  # Simulate !!!
  # Invoke rpoint() for each type separately
  for(i in 1:ntypes) {
    if(verbose) cat(paste("Type", i, "\n"))
    if(single.arg && is.function(f))
      # call f(x,y,m, ...)
      Y <- rpoint(nn[i], f, fmax=maxes[i], win=win,
                  types[i], ..., giveup=giveup, verbose=verbose)
    else
      # call f(x,y, ...) or use other formats
      Y <- rpoint(nn[i], flist[[i]], fmax=maxes[i], win=win,
                  ..., giveup=giveup, verbose=verbose)
    Y <- Y %mark% factortype[i]
    X <- if(i == 1) Y else superimpose(X, Y)
  }
  
  # Randomly permute, in case the order is important
  permu <- sample(X$n)
  return(X[permu])
}

rmpoint.I.allim <- function(n, f, types) {
  # Internal use only!
  # Generates random marked points (Model I *only*)
  # when all f[[i]] are pixel images.
  #
  # Extract pixel coordinates and probabilities
  get.stuff <- function(imag) {
    w <- as.mask(as.owin(imag))
    dx <- w$xstep
    dy <- w$ystep
    xpix <- as.vector(raster.x(w)[w$m])
    ypix <- as.vector(raster.y(w)[w$m])
    ppix <- as.vector(imag$v[w$m]) # not normalised - OK
    npix <- length(xpix)
    return(list(xpix=xpix, ypix=ypix, ppix=ppix,
                dx=rep(dx,npix), dy=rep(dy, npix),
                npix=npix))
  }
  stuff <- lapply(f, get.stuff)
  # Concatenate into loooong vectors
  xpix <- unlist(lapply(stuff, function(z) { z$xpix }))
  ypix <- unlist(lapply(stuff, function(z) { z$ypix }))
  ppix <- unlist(lapply(stuff, function(z) { z$ppix }))
  dx <- unlist(lapply(stuff, function(z) { z$dx }))
  dy <- unlist(lapply(stuff, function(z) { z$dy }))
  # replicate types
  numpix <- unlist(lapply(stuff, function(z) { z$npix }))
  tpix <- rep(seq(types), numpix)
  #
  # sample pixels from union of all images
  #
  npix <- sum(numpix)
  id <- sample(npix, n, replace=TRUE, prob=ppix)
  # get pixel centre coordinates and randomise within pixel
  x <- xpix[id] + (runif(n) - 1/2) * dx[id]
  y <- ypix[id] + (runif(n) - 1/2) * dy[id]
  # compute types
  marx <- factor(types[tpix[id]],levels=types)
  # et voila!
  return(ppp(x, y, window=as.owin(f[[1]]), marks=marx))
}

#
#     wrapper for Rolf's function
#
rpoint.multi <- function (n, f, fmax=NULL, marks = NULL, win =
                    unit.square(), giveup = 1000, verbose = FALSE) {
  # unmarked case
  if (length(marks) <= 1) {
    if(is.function(f))
      return(rpoint(n, f, fmax, win, giveup=giveup, verbose=verbose))
    else
      return(rpoint(n, f, fmax, giveup=giveup, verbose=verbose))
  }
  # multitype case
  if(length(marks) != n)
    stop("length of marks vector != n")
  if(!is.factor(marks))
    stop("marks should be a factor")
  types <- levels(marks)
  types <- factor(types, levels=types)
  # generate required number of points of each type
  nums <- table(marks)
  X <- rmpoint(nums, f, fmax, win=win, types=types,
               giveup=giveup, verbose=verbose)
  if(any(table(X$marks) != nums))
    stop("Internal error: output of rmpoint illegal")
  # reorder them to correspond to the desired 'marks' vector
  Y <- X
  for(ty in types) {
    to   <- (marks == ty)
    from <- (X$marks == ty)
    if(sum(to) != sum(from))
      stop(paste("Internal error: mismatch for mark =", ty))
    if(any(to)) {
      Y$x[to] <- X$x[from]
      Y$y[to] <- X$y[from]
      Y$marks[to] <- ty
    }
  }
  return(Y)
}


  
  

    
