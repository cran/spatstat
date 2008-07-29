#
#    random.S
#
#    Functions for generating random point patterns
#
#    $Revision: 4.36 $   $Date: 2008/07/28 01:13:42 $
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
#    rthin()           independent random thinning
#    rjitter()         random perturbation
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
      return(ppp(x, y, window=win, check=FALSE))
}

"runifdisc" <-
  function(n, radius=1, centre=c(0,0), ...)
{
  # i.i.d. uniform points in the disc of radius r and centre (x,y)
  disque <- disc(centre=centre, radius=radius, ...)
  theta <- runif(n, min=0, max= 2 * pi)
  s <- sqrt(runif(n, min=0, max=radius^2))
  return(ppp(centre[1] + s * cos(theta),
             centre[2] + s * sin(theta),
             window=disque, check=FALSE))
}


"runifpoint" <-
  function(n, win=owin(c(0,1),c(0,1)), giveup=1000)
{
    win <- as.owin(win)
    if(n == 0)
      return(ppp(numeric(0), numeric(0), window=win))
    
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
             return(ppp(x, y, window=win, check=FALSE))
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
               qq <- qq[win]
               # add them to result
               X <- superimpose(X, qq, W=win)
               # if we have enough points, exit
               if(X$n > n) 
                 return(X[1:n])
               else if(X$n == n)
                 return(X)
               # otherwise get bored eventually
               else if(ntries >= giveup)
                 stop(paste("Gave up after", giveup * n, "trials,",
                            X$n, "points accepted"))
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
      X <- X[win]  

    return(X)
}

rpoint <- function(n, f, fmax=NULL,
                   win=unit.square(), ..., giveup=1000,verbose=FALSE) {
  
  if(missing(f) || (is.numeric(f) && length(f) == 1))
    # uniform distribution
    return(runifpoint(n, win, giveup))
  
  # non-uniform distribution....
  
  if(!is.function(f) && !is.im(f))
    stop(paste(sQuote("f"),
               "must be either a function or an",
               sQuote("im"), "object"))
  
  if(is.im(f)) {
    # ------------ PIXEL IMAGE ---------------------
    wf <- as.owin(f)
    if(n == 0)
      return(ppp(numeric(0), numeric(0), window=wf))
    w <- as.mask(wf)
    M <- w$m
    dx <- w$xstep
    dy <- w$ystep
    # extract pixel coordinates and probabilities
    xpix <- as.vector(raster.x(w)[M])
    ypix <- as.vector(raster.y(w)[M])
    ppix <- as.vector(f$v[M]) # not normalised - OK
    # select pixels
    id <- sample(length(xpix), n, replace=TRUE, prob=ppix)
    # extract pixel centres and randomise within pixels
    x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
    y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
    return(ppp(x, y, window=wf, check=FALSE))
  }

  # ------------ FUNCTION  ---------------------  
  # Establish parameters for rejection method

  verifyclass(win, "owin")
  if(n == 0)
    return(ppp(numeric(0), numeric(0), window=win))
  
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
      prop <- prop[win]
    if(prop$n > 0) {
      fvalues <- f(prop$x, prop$y, ...)
      paccept <- fvalues/fmax
      u <- runif(prop$n)
      # accepted points
      Y <- prop[u < paccept]
      if(Y$n > 0) {
        # add to X
        X <- superimpose(X, Y, W=win)
        if(X$n >= n) {
          # we have enough!
          if(verbose)
            cat(paste("acceptance rate = ",
                      round(100 * X$n/(ntries * n), 2), "%\n"))
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

    if(!(is.numeric(lambda) || is.function(lambda) || is.im(lambda)))
      stop(paste(sQuote("lambda"),
                 "must be a constant, a function or an image"))
    if(is.numeric(lambda) && !(length(lambda) == 1 && lambda >= 0))
      stop(paste(sQuote("lambda"),
                 "must be a single, nonnegative number"))
    if(!is.null(lmax)) {
      if(!is.numeric(lmax))
        stop("lmax should be a number")
      if(length(lmax) > 1)
        stop("lmax should be a single number")
    }
      
    win <- if(is.im(lambda))
      rescue.rectangle(as.owin(lambda))
    else
      as.owin(win)
    
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
    stop(paste(sQuote("lambda"), "must be a constant, a function or an image"))
}
    
"rMaternI" <-
  function(kappa, r, win = owin(c(0,1),c(0,1)))
{
    X <- rpoispp(kappa, win=win)
    if(X$n <= 1) return(X)
    d <- nndist(X$x, X$y)
    qq <- X[d > r]
    return(qq)
}
    
"rMaternII" <-
  function(kappa, r, win = owin(c(0,1),c(0,1)))
{
    X <- rpoispp(kappa, win=win)
    
    if(X$n <= 1) return(X)
    
    # matrix of squared pairwise distances
    d2 <- pairdist(X$x, X$y, squared=TRUE)
    close <- (d2 <= r^2)

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
  stopifnot(is.numeric(r) && length(r) == 1 && r >= 0)
  stopifnot(is.numeric(n) && length(n) == 1 && n >= 0)
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
         X <- superimpose(X, qq, W=win)
       if(X$n == n)
         return(X)
     }
     warning(paste("Gave up after", giveup,
                "attempts with only", X$n, "points placed out of", n))
     return(X)
}

"rNeymanScott" <-
  function(kappa, rmax, rcluster, win = owin(c(0,1),c(0,1)), ..., lmax=NULL)
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
  if(is.im(kappa) && !is.subset.owin(as.owin(kappa), dilated))
    stop(paste("The window in which the image",
               sQuote("kappa"),
               "is defined\n",
               "is not large enough to contain the dilation of the window",
               sQuote("win")))
  parents <- rpoispp(kappa, lmax=lmax, win=dilated)
  #
  result <- NULL
  # generate clusters
  if(parents$n > 0) {
    for(i in seq(parents$n)) {
      # generate random offspring of i-th parent point
      cluster <- rcluster(parents$x[i], parents$y[i], ...)
      if(!inherits(cluster, "ppp"))
        cluster <- ppp(cluster$x, cluster$y, window=frame, check=FALSE)
      # skip if cluster is empty
      if(cluster$n > 0) {
        # trim to window
        cluster <- cluster[win]
        if(is.null(result)) {
          # initialise offspring pattern and offspring-to-parent map
          result <- cluster
          parentid <- rep(1, cluster$n)
        } else {
          # add to pattern
          result <- superimpose(result, cluster, W=win)
          # update offspring-to-parent map
          parentid <- c(parentid, rep(i, cluster$n))
        }
      }
    }
  }

  attr(result, "parents") <- parents
  attr(result, "parentid") <- parentid
  
  return(result)
}  

"rMatClust" <-
  function(kappa, r, mu, win = owin(c(0,1),c(0,1)))
{
  # Matern Cluster Process with Poisson (mu) offspring distribution
  #
  stopifnot(is.numeric(r) && length(r) == 1 && r > 0)
  #
  if(is.numeric(mu)) {
    # homogeneous
    stopifnot(length(mu) == 1 && mu >= 0) 
    poisclus <-
      function(x0, y0, radius, mu) {
        n <- rpois(1, mu)
        return(runifdisc(n, centre=c(x0, y0), radius=radius))
      }
    result <- rNeymanScott(kappa, r, poisclus, win, radius=r, mu=mu)
  } else {
    # inhomogeneously modulated clusters a la Waagepetersen
    mu <- as.im(mu, win)
    poisclusinhom <-
      function(x0, y0, radius, mu, mwin, insidedisc) {
        disque <- as.im(function(x,y,x0,y0,rsq) { (x-x0)^2 + (y-y0)^2 < rsq },
                        mwin, x0=x0,y0=y0, rsq=radius^2)
        intens <- eval.im(mu * disque / (pi * radius^2))
        clust <- rpoispp(intens)
      }
    result <- rNeymanScott(kappa, r, poisclusinhom, win, radius=r, mu=mu,
                           mwin=as.owin(mu))
  }
  return(result)
}
    
"rThomas" <-
  function(kappa, sigma, mu, win = owin(c(0,1),c(0,1)))
{
  # Thomas process with Poisson(mu) number of offspring
  # at isotropic Normal(0,sigma^2) displacements from parent
  #
  stopifnot(is.numeric(sigma) && length(sigma) == 1 && sigma > 0)

  if(is.numeric(mu)) {
    # homogeneous
    stopifnot(length(mu) == 1 && mu >= 0)
    thomclus <-  function(x0, y0, sigma, mu) {
                           n <- rpois(1, mu)
                           x <- rnorm(n, mean=x0, sd=sigma)
                           y <- rnorm(n, mean=y0, sd=sigma)
                           return(list(x=x, y=y))
                         }
    result <- rNeymanScott(kappa, 4 * sigma, thomclus,
                         win, sigma=sigma, mu=mu)
  } else {
    # inhomogeneously modulated clusters a la Waagepetersen
    mu <- as.im(mu, win)
    gaussxy <- function(x,y,x0,y0,sigma) 
      dnorm(x, mean=x0, sd=sigma) * dnorm(y, mean=y0, sd=sigma)
    thomclusinhom <-  function(x0, y0, sigma, mu, gaussxy, mwin) {
      gauss <- as.im(gaussxy, mwin, x0=x0,y0=y0,sigma=sigma)
      mumeasure <- eval.im(gauss * mu)
      clust <- rpoispp(mumeasure)
      return(clust)
    }
    result <- rNeymanScott(kappa, 4 * sigma, thomclusinhom,
                           win, sigma=sigma, mu=mu,
                           gaussxy=gaussxy, mwin=as.owin(mu))
  }
  return(result)
}

rGaussPoisson <-
  function(kappa, r, p2, win=owin(c(0,1), c(0,1)))
{
  # Gauss-Poisson process
  oneortwo <- function(x0, y0, radius, p2) {
    if(runif(1) > p2) 
      # one point
      return(list(x=x0, y=y0))
    # two points
    theta <- runif(1, min=0, max=2*pi)
    return(list(x=x0+c(-1,1)*radius*cos(theta),
                y=y0+c(-1,1)*radius*sin(theta)))
  }
  result <- rNeymanScott(kappa, 1.05 * r, oneortwo,
                         win, radius=r/2, p2=p2)
  return(result)
  
}

  
rstrat <- function(win=square(1), nx, ny=nx, k=1) {
  win <- as.owin(win)
  stopifnot(nx >= 1 && ny >= 1)
  stopifnot(k >= 1)
  xy <- stratrand(win, nx, ny, k)
  Xbox <- ppp(xy$x, xy$y, win$xrange, win$yrange, check=FALSE)
  X <- Xbox[win]
  return(X)
}

rsyst <- function(win=square(1), nx, ny=nx, dx=NULL, dy=NULL) {
  win <- as.owin(win)
  xr <- win$xrange
  yr <- win$yrange
  if(!missing(nx)) {
    stopifnot(nx >= 1 && ny >= 1)
    if(!is.null(dx) || !is.null(dy))
      stop("Do not give both nx and dx")
    dx <- diff(xr)/nx
    dy <- diff(yr)/ny
    x0 <- seq(xr[1], xr[2], length=nx+1)
    y0 <- seq(yr[1], yr[2], length=ny+1)
  } else if(!is.null(dx)) {
    stopifnot(dx > 0 & dy > 0)
    if(!is.null(nx) || !is.null(ny))
      stop("Do not give both nx and dx")
    x0 <- seq(xr[1], xr[2], by=dx)
    y0 <- seq(yr[1], yr[2], by=dy)
  } else stop("Either (nx, ny) or (dx, dy) should be given")
  xy0 <- expand.grid(x=x0, y=y0)
  x <- xy0$x + runif(1, min = 0, max = dx)
  y <- xy0$y + runif(1, min = 0, max = dy)
  Xbox <- ppp(x, y, xr, yr, check=FALSE)
  X <- Xbox[win]
  return(X)
}


rcell <- function(win=square(1), nx, ny=nx, dx=NULL, dy=NULL) {
  win <- as.owin(win)
  xr <- win$xrange
  yr <- win$yrange
  if(!missing(nx)) {
    stopifnot(nx >= 1 && ny >= 1)
    if(!is.null(dx) || !is.null(dy))
      stop("Do not give both nx and dx")
    dx <- diff(xr)/nx
    dy <- diff(yr)/ny
    x0 <- seq(xr[1], xr[2], length=nx+1)
    y0 <- seq(yr[1], yr[2], length=ny+1)
  } else if(!is.null(dx)) {
    stopifnot(dx > 0 & dy > 0)
    if(!is.null(nx) || !is.null(ny))
      stop("Do not give both nx and dx")
    x0 <- seq(xr[1], xr[2], by=dx)
    y0 <- seq(yr[1], yr[2], by=dy)
  } else stop("Either (nx, ny) or (dx, dy) should be given")
  mx <- length(x0)
  my <- length(y0)
  x <- numeric(0)
  y <- numeric(0)
  rcellnumber <- function(n) {
    u <- runif(n, min=0, max=1)
    k <- ifelse(u < 1/10, 0, ifelse(u < 89/90, 1, 10))
    return(k)
  }
  for(ix in seq(mx))
    for(iy in seq(my)) {
      nij <- rcellnumber(1)
      x <- c(x, x0[ix] + runif(nij, min=0, max=dx))
      y <- c(y, y0[iy] + runif(nij, min=0, max=dy))
    }
  Xbox <- ppp(x, y, xr, yr, check=FALSE)
  X <- Xbox[win]
  return(X)
}


rthin <- function(X, P, ...) {
  verifyclass(X, "ppp")

  if(is.numeric(P)) {
    # vector of retention probabilities
    pX <- P
    if(length(pX) != X$n) {
      if(length(pX) == 1)
        pX <- rep(pX, X$n)
      else 
        stop("Length of vector P does not match number of points of X")
    }
    if(any(is.na(pX)))
      stop("P contains NA's")
  } else if(is.function(P)) {
    # function - evaluate it at points of X
    pX <- P(X$x, X$y, ...)
    if(length(pX) != X$n)
      stop("Function P returned a vector of incorrect length")
    if(!is.numeric(pX))
      stop("Function P returned non-numeric values")
    if(any(is.na(pX)))
      stop("Function P returned some NA values")
    prange <- range(pX)
  } else if(is.im(P)) {
    # image - look it up
    if(!(P$type %in% c("integer", "real")))
      stop("Values of image P should be numeric")
    pX <- P[X, drop=FALSE]
    if(any(is.na(pX)))
      stop("some points of X lie outside the domain of image P")
  } else
  stop("Unrecognised format for P")

  if(min(pX) < 0) stop("some probabilities are negative")
  if(max(pX) > 1) stop("some probabilities are greater than 1")

  retain <- (runif(length(pX)) < pX)
  return(X[retain])
}


# rjitter

rjitter <- function(X, radius, retry=TRUE, giveup=10000) {
  verifyclass(X, "ppp")
  W <- X$window
  if(!retry) {
    # points outside window are lost
    D <- runifdisc(X$n, radius=radius)
    xnew <- X$x + D$x
    ynew <- X$y + D$y
    ok <- inside.owin(xnew, ynew, W)
    return(ppp(xnew[ok], ynew[ok], window=W))
  }
  # retry = TRUE: condition on points being inside window
  undone <- rep(TRUE, X$n)
  while(any(undone)) {
    giveup <- giveup - 1
    if(giveup <= 0)
      return(X)
    Y <- X[undone]
    D <- runifdisc(Y$n, radius=radius)
    xnew <- Y$x + D$x
    ynew <- Y$y + D$y
    ok <- inside.owin(xnew, ynew, W)
    if(any(ok)) {
      changed <- seq(X$n)[undone][ok]
      X$x[changed] <- xnew[ok]
      X$y[changed] <- ynew[ok]
      undone[ok] <- FALSE
    }
  }
  return(X)
}

