#
#    random.S
#
#    Functions for generating random point patterns
#
#    $Revision: 4.12 $   $Date: 2004/01/27 07:05:16 $
#
#
#    runifpoint()      n i.i.d. uniform random points ("binomial process")
#
#    runifpoispp()     uniform Poisson point process
#
#    rpoispp()         general Poisson point process (thinning method)
#
#    rpoint()          n independent random points (rejection/pixel list)
#
#    rMaternI()        Mat'ern model I 
#    rMaternII()       Mat'ern model II
#    rSSI()            Simple Sequential Inhibition process
#
#    rNeymanScott()    Neyman-Scott process (generic)
#    rMatClust()       Mat'ern cluster process
#    rThomas()         Thomas process
#
#
#
#    Examples:
#          u01 <- owin(0:1,0:1)
#          plot(runifpoispp(100, u01))
#          X <- rpoispp(function(x,y) {100 * (1-x/2)}, 100, u01)
#          X <- rpoispp(function(x,y) {ifelse(x < 0.5, 100, 20)}, 100)
#          plot(X)
#          plot(rMaternI(100, 0.02))
#          plot(rMaternII(100, 0.05))
#

"runifrect" <-
  function(n, win=owin(c(0,1),c(0,1)))
{
  # no checking
      x <- runif(n, min=win$xrange[1], max=win$xrange[2])
      y <- runif(n, min=win$yrange[1], max=win$yrange[2])  
      return(ppp(x, y, window=win))
}

"runifdisc" <-
  function(n, r=1, x=0, y=0)
{
  # i.i.d. uniform points in the disc of radius r and centre (x,y)
  theta <- runif(n, min=0, max= 2 * pi)
  s <- sqrt(runif(n, min=0, max=r^2))
  return(list(x = x + s * cos(theta), y = y + s * sin(theta)))
}


"runifpoint" <-
  function(n, win=owin(c(0,1),c(0,1)), giveup=1000)
{
    win <- as.owin(win)

    switch(win$type,
           rectangle = {
             return(runifrect(n, win))
           },
           mask = {
             dx <- win$xstep
             dy <- win$ystep
             # extract pixel coordinates and probabilities
             xpix <- as.vector(raster.x(win)[win$m])
             ypix <- as.vector(raster.y(win)[win$m])
             # select pixels with equal probability
             id <- sample(seq(xpix), n, replace=TRUE)
             # extract pixel centres and randomise within pixels
             x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
             y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
             return(ppp(x, y, window=win))
           },
           polygonal={
             # rejection method
             # initialise empty pattern
             x <- numeric(0)
             y <- numeric(0)
             X <- ppp(x, y, window=win)
             #
             # rectangle in which trial points will be generated
             box <- bounding.box(win)
             # 
             ntries <- 0
             repeat {
               ntries <- ntries + 1
               # generate trial points in batches of n
               qq <- runifrect(n, box) 
               # retain those which are inside 'win'
               qq <- qq[, win]
               # add them to result
               X <- superimpose(X, qq)
               # if we have enough points, exit
               if(X$n > n) 
                 return(X[1:n])
               else if(X$n == n)
                 return(X)
               # otherwise get bored eventually
               else if(ntries >= giveup)
                 stop(paste("Gave up after", giveup * n, "trials,",
                            np, "points accepted"))
             }
           })
    stop("Unrecognised window type")
}

"runifpoispp" <-
function(lambda, win = owin(c(0,1),c(0,1))) {
    win <- as.owin(win)
    if(!is.numeric(lambda) || length(lambda) > 1 || lambda < 0)
      stop("Intensity lambda must be a single number >= 0")

    if(lambda == 0) # return empty pattern
      return(ppp(numeric(0), numeric(0), window=win))

    # generate Poisson process in enclosing rectangle 
    box <- bounding.box(win)
    mean <- lambda * area.owin(box)
    n <- rpois(1, mean)
    X <- runifpoint(n, box)

    # trim to window
    if(win$type != "rectangle")
      X <- X[, win]  

    return(X)
}

rpoint <- function(n, f, fmax=NULL,
                   win=unit.square(), ..., giveup=1000,verbose=FALSE) {
  
  if(missing(f) || (is.numeric(f) && length(f) == 1))
    # uniform distribution
    return(runifpoint(n, win, giveup))
  
  # non-uniform distribution....
  
  if(!is.function(f) && !is.im(f))
    stop("\`f\' must be either a function or an \`im\' object")
  
  if(is.im(f)) {
    # ------------ PIXEL IMAGE ---------------------
    w <- as.mask(as.owin(f))
    dx <- w$xstep
    dy <- w$ystep
    # extract pixel coordinates and probabilities
    xpix <- as.vector(raster.x(w)[w$m])
    ypix <- as.vector(raster.y(w)[w$m])
    ppix <- as.vector(f$v[w$m]) # not normalised - OK
    # select pixels
    id <- sample(length(xpix), n, replace=TRUE, prob=ppix)
    # extract pixel centres and randomise within pixels
    x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
    y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
    return(ppp(x, y, window=as.owin(f)))
  }

  # ------------ FUNCTION  ---------------------  
  # Establish parameters for rejection method

  verifyclass(win, "owin")
  if(is.null(fmax)) {
    # compute approx maximum value of f
    imag <- as.im(f, win, ...)
    summ <- summary(imag)
    fmax <- summ$max + 0.05 * diff(summ$range)
  }
  irregular <- (win$type != "rectangle")
  box <- bounding.box(win)
  X <- ppp(numeric(0), numeric(0), window=win)
  
  ntries <- 0

  # generate uniform random points in batches
  # and apply the rejection method.
  # Collect any points that are retained in X

  repeat{
    ntries <- ntries + 1
    # proposal points
    prop <- runifrect(n, box)
    if(irregular)
      prop <- prop[, win]
    if(prop$n > 0) {
      fvalues <- f(prop$x, prop$y, ...)
      paccept <- fvalues/fmax
      u <- runif(prop$n)
      # accepted points
      Y <- prop[u < paccept]
      if(Y$n > 0) {
        # add to X
        X <- superimpose(X, Y)
        if(X$n >= n) {
          # we have enough!
          if(verbose)
            cat(paste("acceptance rate = ",
                      round(100 * X$n/(ntries * n), 2), "\%\n"))
          return(X[1:n])
        }
      }
    }
    if(ntries > giveup)
      stop(paste("Gave up after",giveup * n,"trials with",
                 X$n, "points accepted"))
  }
  invisible(NULL)
}

"rpoispp" <-
  function(lambda, lmax=NULL, win = owin(c(0,1),c(0,1)), ...) {
    # arguments:
    #     lambda  intensity: constant, function(x,y,...) or image
    #     lmax     maximum possible value of lambda(x,y,...)
    #     win     default observation window (of class 'owin')
    #   ...       arguments passed to lambda(x, y, ...)

    win <- if(is.im(lambda)) as.owin(lambda) else as.owin(win)
    
    if(is.numeric(lambda)) 
      # uniform Poisson
      return(runifpoispp(lambda, win))

    # inhomogeneous Poisson
    # perform thinning of uniform Poisson

    if(is.null(lmax)) {
      imag <- as.im(lambda, win, ...)
      summ <- summary(imag)
      lmax <- summ$max + 0.05 * diff(summ$range)
    }

    if(is.function(lambda)) {
      X <- runifpoispp(lmax, win)  # includes sanity checks on `lmax'
      if(X$n == 0) return(X)
      prob <- lambda(X$x, X$y, ...)/lmax
      u <- runif(X$n)
      retain <- (u <= prob)
      X <- X[retain, ]
      return(X)
    }
    if(is.im(lambda)) {
      X <- runifpoispp(lmax, win)
      if(X$n == 0) return(X)
      prob <- lambda[X]/lmax
      u <- runif(X$n)
      retain <- (u <= prob)
      X <- X[retain, ]
      return(X)
    }
    stop("\'lambda\' must be a constant, a function or an image")
}
    
"rMaternI" <-
  function(lambda, r, win = owin(c(0,1),c(0,1)))
{
    X <- rpoispp(lambda, win=win)
    if(X$n <= 1) return(X)
    d <- nndist(X$x, X$y)
    qq <- X[d > r]
    return(qq)
}
    
"rMaternII" <-
  function(lambda, r, win = owin(c(0,1),c(0,1)))
{
    X <- rpoispp(lambda, win=win)
    
    if(X$n <= 1) return(X)
    
    # matrix of pairwise distances
    d <- pairdist(X$x, X$y)
    close <- (d <= r)

    # random order 1:n
    age <- sample(seq(X$n), X$n, replace=FALSE)
    earlier <- outer(age, age, ">")

    conflict <- close & earlier
    # delete <- apply(conflict, 1, any)
    delete <- matrowany(conflict)
    
    qq <- X[ !delete]
    return(qq)
}
  
"rSSI" <-
  function(r, n, win = owin(c(0,1),c(0,1)), giveup = 1000)
{
     # Simple Sequential Inhibition process
     # fixed number of points
     # Naive implementation, proposals are uniform
     win <- as.owin(win)
     X <- ppp(numeric(0),numeric(0), window=win)
     r2 <- r^2
     if(n * pi * r2/4  > area.owin(win))
       stop(paste("Window is too small to fit", n, "points",
                  "at minimum separation", r))
     ntries <- 0
     while(ntries < giveup) {
       ntries <- ntries + 1
       qq <- runifpoint(1, win)
       x <- qq$x[1]
       y <- qq$y[1]
       if(X$n == 0 || all(((x - X$x)^2 + (y - X$y)^2) > r2))
         X <- superimpose(X, qq)
       if(X$n == n)
         return(X)
     }
     warning(paste("Gave up after", giveup,
                "attempts with only", X$n, "points placed out of", n))
     return(X)
}

"rNeymanScott" <-
  function(lambda, rmax, rcluster, win = owin(c(0,1),c(0,1)), ..., lmax=NULL)
{
  # Generic Neyman-Scott process
  # Implementation for bounded cluster radius
  #
  # 'rcluster' is a function(x,y) that takes the coordinates
  # (x,y) of the parent point and generates a list(x,y) of offspring
  #
  # "..." are arguments to be passed to 'rcluster()'
  #

  win <- as.owin(win)
  
  # Generate parents in dilated window
  frame <- bounding.box(win)
  dilated <- owin(frame$xrange + c(-rmax, rmax),
                  frame$yrange + c(-rmax, rmax))
  if(is.im(lambda) && !is.subset.owin(as.owin(lambda), dilated))
    stop("The window in which the image \'lambda\' is defined\n\
is not large enough to contain the dilation of the window \'win\'")
  parents <- rpoispp(lambda, lmax=lmax, win=dilated)
  #
  result <- ppp(numeric(0), numeric(0), window = win)
  
  if(parents$n == 0)
    return(result)
  
  for(i in seq(parents$n)) {
    # generate random offspring of i-th parent point
    cluster <- rcluster(parents$x[i], parents$y[i], ...)
    cluster <- ppp(cluster$x, cluster$y, window=frame)
    # trim to window
    cluster <- cluster[,win]
    # add to pattern
    result <- superimpose(result, cluster)
  }
  return(result)
}  

"rMatClust" <-
  function(lambda, r, mu, win = owin(c(0,1),c(0,1)))
{
  # Matern Cluster Process with Poisson (mu) offspring distribution
  #
  poisclus <-  function(x0, y0, radius, mu) {
                           n <- rpois(1, mu)
                           return(runifdisc(n, radius, x0, y0))
                         }
  result <- rNeymanScott(lambda, r, poisclus,
                         win, radius=r, mu=mu)
  return(result)
}
    
"rThomas" <-
  function(lambda, sigma, mu, win = owin(c(0,1),c(0,1)))
{
  # Thomas process with Poisson(mu) number of offspring
  # at isotropic Normal(0,sigma^2) displacements from parent
  #
  thomclus <-  function(x0, y0, sigma, mu) {
                           n <- rpois(1, mu)
                           x <- rnorm(n, mean=x0, sd=sigma)
                           y <- rnorm(n, mean=y0, sd=sigma)
                           return(list(x=x, y=y))
                         }
  result <- rNeymanScott(lambda, 4 * sigma, thomclus,
                         win, sigma=sigma, mu=mu)
  return(result)
}
  
