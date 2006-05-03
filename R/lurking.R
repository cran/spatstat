# Lurking variable plot for arbitrary covariate.
#
#
# $Revision: 1.4 $ $Date: 2006/04/24 11:24:10 $
#

lurking <- function(object, covariate, type="eem",
                    cumulative=TRUE,
                    clipwindow=default.clipwindow(object),
                    rv = NULL,
                    plot.sd, plot.it=TRUE,
                    typename,
                    covname, ...) {
  # validate
  verifyclass(object, "ppm")

  # extract spatial locations 
  Q <- quad.ppm(object)
  datapoints <- Q$data
  quadpoints <- union.quad(Q)
  Z <- is.data(Q)
  wts <- w.quad(Q)

  #################################################################
  # compute the covariate

  if(is.vector(covariate) && is.numeric(covariate)) 
    covvalues <- covariate
  else if(is.im(covariate))
    covvalues <- covariate[quadpoints]
  else if(is.expression(covariate)) {
    # Evaluate the desired covariate at all quadrature points
    # Set up environment for evaluating expression
    glmdata <- object$internal$glmdata
    # Fix special cases
    if(is.null(glmdata)) {
      # default 
      glmdata <- data.frame(x=quadpoints$x, y=quadpoints$y)
      if(is.marked(quadpoints))
        glmdata$marks <- quadpoints$marks
    }
    # ensure x and y are in data frame 
    if(!all(c("x","y") %in% names(glmdata))) {
      glmdata$x <- quadpoints$x
      glmdata$y <- quadpoints$y
    }
    # Evaluate expression
    evaluate <- function(a, b) {
      if (exists("is.R") && is.R())
        eval(a, envir = b)
      else eval(a, local = b)
    }
    covvalues <- evaluate(covariate, glmdata)
    if(!is.numeric(covvalues))
      stop("The evaluated covariate is not numeric")
  } else 
    stop("The \`covariate\' should be either a numeric vector or an expression")

  # Validate covariate values
  if(any(is.na(covvalues)))
    stop("covariate contains NA's")
  if(any(is.infinite(covvalues) | is.nan(covvalues)))
    stop("covariate contains Inf or NaN values")
  if(length(covvalues) != quadpoints$n)
    stop("Length of covariate vector,", length(covvalues), "!=",
         quadpoints$n, ", number of quadrature points")

  # Quadrature points marked by covariate value
  covq <- quadpoints %mark% covvalues

  ################################################################
  # Residuals/marks attached to appropriate locations.
  # Stoyan-Grabarnik weights are attached to the data points only.
  # Others (residuals) are attached to all quadrature points.

  resvalues <- 
    if(!is.null(rv)) rv
    else if(type=="eem") eem(object)
    else residuals.ppm(object,type=type)
  
  res <- (if(type == "eem") datapoints else quadpoints) %mark% resvalues

  # ... and the same locations marked by the covariate
  covres <- if(type == "eem") covq[Z] else covq

  # NAMES OF THINGS
  # name of the covariate
  if(missing(covname)) 
    covname <- if(is.expression(covariate)) paste(covariate) else "covariate"
  # type of residual/mark
  if(missing(typename)) 
    typename <- if(!is.null(rv)) "rv" else attr(resvalues, "typename")

  #######################################################################
  # START ANALYSIS
  # Clip to subwindow if needed
  clip <- !is.poisson.ppm(object) ||
              (!missing(clipwindow) && !is.null(clipwindow))
  if(clip) {
    covq <- covq[, clipwindow]
    res <- res[, clipwindow]
    covres <- covres[, clipwindow]
    clipquad <- inside.owin(quadpoints$x, quadpoints$y, clipwindow)
    wts <- wts[ clipquad ]
  }

  # -----------------------------------------------------------------------
  # (A) EMPIRICAL CUMULATIVE FUNCTION
  # based on data points if type="eem", otherwise on quadrature points

    # cumulative sums which ignore NA's
    cumsumna <- function(x) {
      x[is.na(x)] <- 0
      return(cumsum(x))
    }

      # Reorder the data/quad points in order of increasing covariate value
      # and then compute the cumulative sum of their residuals/marks
    o <- order(covres$marks)
    covsort <- covres$marks[o]
    cummark <- cumsumna(res$marks[o])
      # we'll plot(covsort, cummark) in the cumulative case

  # (B) THEORETICAL MEAN CUMULATIVE FUNCTION
  # based on all quadrature points
    
      # Range of covariate values
    covrange <- range(covq$marks, na.rm=TRUE)
      # Suitable breakpoints
    cvalues <- seq(covrange[1], covrange[2], length=100)
    csmall <- cvalues[1] - diff(cvalues[1:2])
    cbreaks <- c(csmall, cvalues)
      # cumulative area as function of covariate values
    covclass <- cut(covq$marks, breaks=cbreaks)
    increm <- tapply(wts, covclass, sum)
    cumarea <- cumsumna(increm)
      # compute theoretical mean (when model is true)
    mean0 <- if(type == "eem") cumarea else rep(0, length(cumarea))
      # we'll plot(cvalues, mean0) in the cumulative case

  # (A'),(B') DERIVATIVES OF (A) AND (B)
  #  Required if cumulative=FALSE  
  #  Estimated by spline smoothing
    if(!cumulative) {
      # fit smoothing spline to (A) 
      ss <- smooth.spline(covsort, cummark, ...)
      # estimate derivative of (A)
      derivmark <- predict(ss, covsort, deriv=1)$y 
      # similarly for (B) 
      ss <- smooth.spline(cvalues, mean0, ...)
      derivmean <- predict(ss, cvalues, deriv=1)$y
    }
  
  # -----------------------------------------------------------------------
  # Store what will be plotted
  
   if(cumulative) {
     empirical <- data.frame(covariate=covsort, value=cummark)
     theoretical <- data.frame(covariate=cvalues, mean=mean0)
   } else {
     empirical <- data.frame(covariate=covsort, value=derivmark)
     theoretical <- data.frame(covariate=cvalues, mean=derivmean)
   }

  # ------------------------------------------------------------------------
  
    # (C) STANDARD DEVIATION if desired
    # (currently implemented only for Poisson)
    # (currently implemented only for cumulative case)

    if(missing(plot.sd))
      plot.sd <- is.poisson.ppm(object)
    else if(plot.sd && !is.poisson.ppm(object))
      warning("standard deviation is calculated for Poisson model; not valid for this model")

    if(plot.sd && cumulative)
      switch(type,
             pearson={
               theoretical$sd <- sqrt(cumarea)
             },
             raw={
               # Compute fitted conditional intensity at quadrature points
               lambda <- fitted.ppm(object, type="trend")
               if(clip) lambda <- lambda[clipquad]
               # Compute sum of w*lambda for quadrature points in each interval
               dvar <- tapply(wts * lambda, covclass, sum)
               # tapply() returns NA when the table is empty
               dvar[is.na(dvar)] <- 0
               # Cumulate
               theoretical$sd <- sqrt(cumsum(dvar))
             },
             inverse=, # same as eem
             eem={
               # Compute fitted conditional intensity at quadrature points
               lambda <- fitted.ppm(object, type="trend")
               if(clip) lambda <- lambda[clipquad]
               # Compute sum of w/lambda for quadrature points in each interval
               dvar <- tapply(wts / lambda, covclass, sum)
               # tapply() returns NA when the table is empty
               dvar[is.na(dvar)] <- 0
               # Cumulate
               theoretical$sd <- sqrt(cumsum(dvar))
             }
     )

    # ---------------  PLOT THEM  ----------------------------------
    if(plot.it) {
      # work out plot range
      mr <- range(c(0, empirical$value, theoretical$mean), na.rm=TRUE)
      if(!is.null(theoretical$sd))
        mr <- range(c(mr, theoretical$mean + 2 * theoretical$sd,
                          theoretical$mean - 2 * theoretical$sd),
                    na.rm=TRUE)

      # start plot
      plot(covrange, mr, type="n", xlab=covname,
           ylab=paste(if(cumulative)"cumulative" else "marginal", typename))
      # (A)/(A') Empirical
      lines(value ~ covariate, empirical)
      # (B)/(B') Theoretical mean
      lines(mean ~ covariate, theoretical, lty=2)
      # (C) Standard deviation 
      if(!is.null(theoretical$sd)) {
        lines(mean + 2 * sd ~ covariate, theoretical, lty=3)
        lines(mean - 2 * sd ~ covariate, theoretical, lty=3)
      }
    }
  
    # ----------------  RETURN COORDINATES ----------------------------
  stuff <- list(empirical=empirical, theoretical=theoretical)

  return(invisible(stuff))
}
