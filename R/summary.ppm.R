#
#    summary.ppm.R
#
#    summary() method for class "ppm"
#
#    $Revision: 1.23 $   $Date: 2007/01/12 01:29:29 $
#
#    summary.ppm()
#    print.summary.ppm()
#
summary.ppm <- function(object, ..., quick=FALSE) {
  verifyclass(object, "ppm")

  x <- object
  y <- list()

  #######  Extract main data components #########################

  QUAD <- object$Q
  DATA <- QUAD$data
  TREND <- x$trend

  INTERACT <- x$interaction
  if(is.null(INTERACT)) INTERACT <- Poisson()

  #######  Check version #########################

  mpl.ver <- versionstring.ppm(object)
  int.ver <- versionstring.interact(INTERACT)
  current <- versionstring.spatstat()

  virgin <- min(package_version(c(mpl.ver, int.ver)))

  y$antiquated <- antiquated <- (virgin <= package_version("1.5"))
  y$old        <- old        <- (virgin < majorminorversion(current))

  y$version    <- as.character(virgin)
                
  ####### Determine type of model ############################
  
  y$no.trend <- identical.formulae(TREND, NULL) || identical.formulae(TREND, ~1)

  y$stationary <- y$no.trend || identical.formulae(TREND, ~marks)

  y$poisson <- is.poisson.interact(INTERACT)

  y$marked <- is.marked.ppp(DATA)
  y$multitype <- is.multitype.ppp(DATA)
  if(y$marked) y$entries <- list(marks = marks(DATA))

  y$name <- paste(
          if(y$stationary) "Stationary " else "Nonstationary ",
          if(y$poisson) {
            if(y$multitype) "multitype "
            else if(y$marked) "marked "
            else ""
          },
          INTERACT$name,
          sep="")

  y$method <- x$method

  y$problems <- x$problems
  
  if(is.logical(quick) && quick) {
    class(y) <- "summary.ppm"
    return(y)
  }
  
  ######  Does it have external covariates?  ####################

  if(!antiquated) {
    covars <- x$covariates
    hc <- !is.null(covars) && (length(covars) > 0)
  } else {
    # Antiquated format
    # Interpret the function call instead
    callexpr <- parse(text=x$call)
    callargs <- names(as.list(callexpr[[1]]))
    # Data frame of covariates was called 'data' in versions up to 1.4-x
    hc <- !is.null(callargs) && !is.na(pmatch("data", callargs))
  }
  y$has.covars <- hc
    
  ######  Arguments in call ####################################
  
  y$args <- x[c("call", "correction", "rbord")]
  
  #######  Main data components #########################

  y$entries <- append(list(quad=QUAD,
                           data=DATA,
                           interaction=INTERACT),
                      y$entries)

  ####### Summarise data ############################

  y$data <- summary.ppp(DATA, checkdup=FALSE)
  y$quad <- summary.quad(QUAD, checkdup=FALSE)

  if(is.character(quick) && (quick == "no prediction"))
    return(y)
  
  ######  Extract fitted model coefficients #########################

  if(exists("is.R") && is.R()) 
    theta <- x$coef # result of coef(glm(...))
  else
    theta <- x$theta # result of dummy.coef(glm(....))

  y$entries$theta <- theta

  # corresponding internal names of regressor variables 
  y$entries$Vnames <- Vnames <- x$internal$Vnames

  ######  Trend component #########################

  y$trend <- list()

  y$trend$name <- if(y$poisson) "Intensity" else "Trend"

  y$trend$formula <- if(y$no.trend) NULL else TREND

  if(y$poisson && y$no.trend) {
    lambda <- exp(theta[[1]])
    if(!y$marked) { 
      y$trend$label <- "Uniform intensity"
      y$trend$value <- lambda
    } else {
      y$trend$label <- "Uniform intensity for each mark level"
      y$trend$value <- lambda
    }
  } else # process is at least one of: marked, nonstationary, non-poisson
  if(y$stationary) {
    if(!y$marked) {
      # stationary non-poisson non-marked
      y$trend$label <- "First order term"
      y$trend$value <- c(beta=exp(theta[[1]]))
    } else {
      # stationary, marked
      mrk <- marks(DATA)
      y$trend$label <-
        if(y$poisson) "Fitted intensities" else "Fitted first order terms"
      # Use predict.ppm to evaluate the fitted intensities
      lev <- factor(levels(mrk), levels=levels(mrk))
      nlev <- length(lev)
      marx <- list(x=rep(0, nlev), y=rep(0, nlev), marks=lev)
      betas <- predict(x, locations=marx, type="trend")
      names(betas) <- paste("beta_", as.character(lev), sep="")
      y$trend$value <- betas
    }
  } else {
    # not stationary 
    y$trend$label <- "Fitted coefficients for trend formula"
    # extract trend terms without trying to understand them much
    if(is.null(Vnames)) 
      trendbits <- theta
    else {
      agree <- outer(names(theta), Vnames, "==")
      whichbits <- apply(!agree, 1, all)
      trendbits <- theta[whichbits]
    }
    # decide whether there are 'labels within labels'
    unlabelled <- unlist(lapply(trendbits,
                                function(x) { is.null(names(x)) } ))
    if(all(unlabelled))
      y$trend$value <- unlist(trendbits)
    else {
      y$trend$value <- list()
      for(i in seq(trendbits))
          y$trend$value[[i]] <-
            if(unlabelled[i])
              unlist(trendbits[i])
            else
              trendbits[[i]]
    }
  }
  
  ######  Interaction component #########################

  if(!y$poisson) {
    if(!is.null(INTERACT$interpret)) {
      # invoke auto-interpretation feature
      if(newstyle.coeff.handling(INTERACT))
        sensible <- (INTERACT$interpret)(x$coef[Vnames], INTERACT)
      else 
        sensible <- (INTERACT$interpret)(x$coef, INTERACT)
      header <- paste("Fitted", sensible$inames)
      printable <- sensible$printable
    } else {
      # fallback
      sensible <- NULL
      header <- "Fitted interaction terms"
      printable <-  exp(unlist(theta[Vnames]))
    }
    y$interaction <- list(sensible=sensible,
                          header=header,
                          printable=printable)
  }

  class(y) <- "summary.ppm"
  return(y)
}

print.summary.ppm <- function(x, ...) {

  if(x$old)
    warning("Model was fitted by an older version of spatstat")
  
  if(is.null(x$args)) {
    # this is the quick version
    cat(paste(x$name, "\n"))
    return(invisible(NULL))
  }

  # otherwise - full details
  cat("Point process model\n")
  howfitted <-
    if(is.null(x$method))
      "unspecified method"
    else
      switch(x$method,
             mpl="maximum pseudolikelihood (Berman-Turner approximation)",
             ho="Huang-Ogata method (approximate maximum likelihood)",
             paste("unrecognised method", sQuote(x$method)))
  cat(paste("fitted by", howfitted, "\n"))

  cat("Call:\n")
  print(x$args$call)

  if(x$old) 
    cat(paste("** Executed by old spatstat version", x$version, " **\n"))
  
  cat(paste("Edge correction:", dQuote(x$args$correction)))
  if(x$args$rbord > 0)
    cat(paste("border correction distance r =", x$args$rbord,"\n"))

  cat("\n----------------------------------------------------\n")

  # print summary of quadrature scheme
  print(x$quad)
  
  cat("\n----------------------------------------------------\n")
  cat("FITTED MODEL:\n\n")

  # This bit is currently identical to print.ppm()
  # except for a bit more fanfare
  # and the inclusion of the 'gory details' bit
  
  notrend <-    x$no.trend
  stationary <- x$stationary
  poisson <-    x$poisson
  markeddata <- x$marked
  multitype  <- x$multitype
        
  markedpoisson <- poisson && markeddata

  # names of interaction variables if any
  Vnames <- x$entries$Vnames
  # their fitted coefficients
  theta <- x$theta

  # ----------- Print model type -------------------
        
  cat(x$name)
  cat("\n")

  if(markeddata) mrk <- x$entries$marks
  if(multitype) {
    cat("Possible marks: \n")
    cat(paste(levels(mrk)))
  }

  # ----- trend --------------------------

  cat(paste("\n\n ---- ", x$trend$name, ": ----\n\n", sep=""))

  if(!notrend) {
    cat("Trend formula: ")
    print(x$trend$formula)
    if(x$has.covars)
      cat("Model involves external covariates\n")
  }
        
  cat(paste("\n", x$trend$label, ":\n", sep=""))
  
  tv <- x$trend$value
  if(!is.list(tv))
    print(tv)
  else 
    for(i in seq(tv))
      print(tv[[i]])
        
  # ---- Interaction ----------------------------

  if(!poisson) {
    cat("\n\n ---- Interaction: -----\n\n")
    print(x$entries$interaction)
    
    cat(paste(x$interaction$header, ":\n", sep=""))
    print(x$interaction$printable)
  }

  ####### Gory details ###################################
  cat("\n\n----------- gory details -----\n")
  theta <- x$entries$theta
      
  cat("\nFitted regular parameters (theta): \n")
  print(theta)

  cat("\nFitted exp(theta): \n")
  print(exp(unlist(theta)))

  ##### Warnings issued #######

  probs <- x$problems
  if(!is.null(probs) && is.list(probs) && (length(probs) > 0)) 
    lapply(probs,
           function(a) {
             if(is.list(a) && !is.null(p <- a$print))
               cat(paste("Problem:\n", p, "\n\n"))
           })
          
  return(invisible(NULL))
}

no.trend.ppm <- function(x) {
  summary.ppm(x, quick=TRUE)$no.trend
}

is.stationary.ppm <- function(x) {
  summary.ppm(x, quick=TRUE)$stationary
}

is.poisson.ppm <- function(x) {
  summary.ppm(x, quick=TRUE)$poisson
}

is.marked.ppm <- function(X, ...) {
  summary.ppm(X, quick=TRUE)$marked
}

is.multitype.ppm <- function(X, ...) {
  summary.ppm(X, quick=TRUE)$multitype
}

