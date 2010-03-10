#
#  effectfun.R
#
#   $Revision: 1.5 $ $Date: 2009/07/16 00:37:02 $
#

effectfun <- function(model, covname, ...) {
  stopifnot(is.ppm(model))
  # determine names of covariates involved
  intern.names <-
    if(is.marked.ppm(model)) c("x", "y", "marks") else c("x", "y")
  co <- model$covariates
  extern.names <- names(co)
  # find the relevant covariate 
  if(missing(covname))
    stop("covname must be provided")
  if(!(covname %in% c(intern.names, extern.names)))
    stop(paste("model does not have a covariate called", sQuote(covname)))
  # check that fixed values for all other covariates are provided 
  given.covs <- names(list(...))
  if(any(uhoh <- !(extern.names %in% c(given.covs, covname)))) {
    nuh <- sum(uhoh)
    stop(paste(ngettext(nuh,
                        "A value for the covariate",
                        "Values for the covariates"),
               commasep(dQuote(extern.names[uhoh])),
               "must be provided (as",
               ngettext(nuh, "an argument", "arguments"),
               "to effectfun)"))
  }
  # establish type and range of covariate values
  N0 <- 256
  if(covname == "x") {
    covtype <- "real"
    W <- as.owin(data.ppm(model))
    Zr <- W$xrange
    Zvals <- seq(Zr[1], Zr[2], length=N0)
  } else if(covname == "y") {
    covtype <- "real"
    W <- as.owin(data.ppm(model))
    Zr <- W$yrange
    Zvals <- seq(Zr[1], Zr[2], length=N0)
  } else if(covname == "marks") {
    covtype <- "factor"
    Zvals <- levels(marks(data.ppm(model)))
  } else {
    # covariate is external
    if(is.data.frame(co)) {
      Z <- co$covname
      covtype <- typeof(Z)
      if(covtype == "double")
        covtype <- "real"
      switch(covtype,
             real={
               Zr <- range(Z)
               Zvals <- seq(Zr[1], Zr[2], length=N0)
             },
             integer={
               Zr <- range(Z)
               Zvals <- seq(Zr[1], Zr[2], by=ceiling((diff(Zr)+1)/N0))
             },
             factor={
               Zvals <- levels(Z)
             },
             logical={
               Zvals <- c(FALSE, TRUE)
             },
             stop(paste("Cannot handle covariate of type", dQuote(covtype)))
             )
    } else if(is.list(co) && all(unlist(lapply(co, is.im)))) {
      Z <- co[[covname]]
      covtype <- Z$type
      switch(covtype,
             real={
               Zr <- summary(Z)$range
               Zvals <- seq(Zr[1], Zr[2], length=N0)
             },
             factor={
               Zvals <- levels(Z)
             },
             logical={
               Zvals <- c(FALSE, TRUE)
             },
             stop(paste("Cannot handle covariate of type", dQuote(covtype)))
             )
    } else stop("Unrecognised format for covariates in model")
  }
  # set up data frames of fake data for predict method
  # First set up default, constant value for each covariate
  N <- length(Zvals)
  fakeloc <- resolve.defaults(list(...),
                              list(x=0, y=0))[c("x","y")]
  if(is.marked.ppm(model)) {
    lev <- levels(marks(data.ppm(model)))
    fakeloc$marks <- lev[1]
  }
  fakeloc <- lapply(fakeloc, function(x,N) { rep(x[1],N)}, N=N)
  fakecov <- lapply(list(...), function(x,N) { rep(x[1],N)}, N=N)
  # Overwrite value for covariate of interest
  if(covname %in% intern.names)
    fakeloc[[covname]] <- Zvals
  else fakecov[[covname]] <- Zvals
  # convert to data frame
  fakeloc <- do.call("data.frame", fakeloc)
  fakecov <- do.call("data.frame", fakecov)
  #
  # Now predict
  lambda <- predict(model, locations=fakeloc, covariates=fakecov)
  #
  dfin <- cbind(fakeloc, fakecov)[covname]
  df <- cbind(dfin, data.frame(lambda=lambda))
  if(covtype == "real") {
    f <- fv(df, argu=covname, ylab="lambda", valu="lambda", alim=Zr,
            desc=c(paste("value of covariate", covname),
              "fitted intensity"))
    return(f)
  } else return(df)
}
