#
# profilepl.R
#
#  $Revision: 1.8 $  $Date: 2010/05/07 11:08:25 $
#
#  computes profile log pseudolikelihood
#

profilepl <- function(s, f, ..., rbord=NULL, verbose=TRUE) {
  s <- as.data.frame(s)
  n <- nrow(s)
  fname <- deparse(substitute(f))
  stopifnot(is.function(f))
  fargs <- names(formals(f))
  parms <- names(s)
  if(!all(fargs %in% parms))
    stop("Some arguments of f are not provided in s")
  if(!all(parms %in% fargs))
    stop("Some variables in s are not arguments of f")
  interlist <- list()
  logmpl <- numeric(n)
  # make a fake call
  pseudocall <- match.call()
  pseudocall[[1]] <- as.symbol("ppm")
  namcal <- names(pseudocall)
  # remove 's' argument
  retain <- (namcal != "s")
  pseudocall <- pseudocall[retain]
  namcal <- namcal[retain]
  # place 'f' argument third 
  np <- length(pseudocall)
  fpos <- (1:np)[namcal == "f"]
  indices <- (1:np)[-fpos]
  if(length(indices) < 3)
    indices <- c(indices, fpos)
  else 
    indices <- c(indices[1:3], fpos, indices[-(1:3)])
  pseudocall <- pseudocall[indices]
  namcal <- names(pseudocall)
  namcal[namcal=="f"] <- "interaction"
  names(pseudocall) <- namcal
  # determine border correction distance
  if(missing(rbord)) {
    # compute rbord = max reach of interactions
    if(verbose) message("(computing rbord)")
    for(i in 1:n) {
      fi <- do.call("f", as.list(s[i,]))
      if(!inherits(fi, "interact"))
        stop(paste("f did not yield an object of class", sQuote("interact")))
      re <- reach(fi)
      if(is.null(rbord))
        rbord <- re
      else if(rbord < re)
        rbord <- re
    }
  } 
  if(verbose) message(paste("fitting", n, "models..."))
  for(i in 1:n) {
    if(verbose)
      progressreport(i, n)
    fi <- do.call("f", as.list(s[i,]))
    if(!inherits(fi, "interact"))
      stop(paste("f did not yield an object of class", sQuote("interact")))
    # fit model
    if(i == 1) {
      # fit from scratch
      fiti <- ppm(interaction=fi, ..., rbord=rbord, savecomputed=TRUE)
      # save intermediate computations (pairwise distances, etc)
      precomp <- fiti$internal$computed
    } else {
      # use precomputed data
      fiti <- ppm(interaction=fi, ..., rbord=rbord, precomputed=precomp)
    }
    # save log PL for each fit
    logmpl[i] <- as.numeric(logLik(fiti))
    # save fitted coefficients for each fit
    co <- coef(fiti)
    if(i == 1) {
      allcoef <- data.frame(matrix(co, nrow=1))
      names(allcoef) <- names(co)
    } else
      allcoef <- rbind(allcoef, co)
  }
  if(verbose) message("done.")
  opti <- which.max(logmpl)
  optint <- do.call("f", as.list(s[opti,, drop=FALSE]))
  optfit <- do.call("ppm",
                    list(interaction=optint, ..., rbord=rbord))
  result <- list(param=s,
                 prof=logmpl,
                 iopt=opti,
                 fit=optfit,
                 rbord=rbord,
                 fname=fname,
                 allcoef=allcoef,
                 otherstuff=list(...),
                 pseudocall=pseudocall)
  class(result) <- c("profilepl", class(result))
  return(result)
}

#
#   print method
#

print.profilepl <- function(x, ...) {
  cat("Profile log pseudolikelihood values\n")
  cat("for model:\t")
  print(x$pseudocall)
  cat(paste("fitted with rbord=", x$rbord, "\n"))
  nparm <- ncol(x$param)
  cat(paste("Interaction:", x$fname,
            "\n", "with irregular",
            ngettext(nparm, "parameter ", "parameters\n")))
  for(na in names(x$param)) {
    ra <- range(x$param[[na]])
    cat(paste(sQuote(na), "in",
              paste("[",ra[1],", ",ra[2],"]",sep=""),
              "\n"))
  }
  cat(paste("Optimum",
            ngettext(nparm, "value", "values"),
            "of irregular",
            ngettext(nparm, "parameter: ", "parameters:\n")))
  popt <- x$param[x$iopt,, drop=FALSE]
  cat(commasep(paste(names(popt), "=", popt)))
  cat("\n")
}

#
#   summary method
#

summary.profilepl <- function(object, ...) {
  print(object)
  cat("\n\nOptimal model:\n")
  print(object$fit)
}


#
#  plot method 
#

plot.profilepl <- function(x, ..., add=FALSE, main=NULL,
                           tag=TRUE, coeff=NULL, xvariable=NULL) {
  para <- x$param
  npara <- ncol(para)
  # main header
  if(is.null(main))
    main <- deparse(x$pseudocall)
  # x variable for plot
  if(is.null(xvariable)) {
    xvalues <- para[,1]
    xname <- names(para)[1]
  } else {
    stopifnot(is.character(xvariable))
    if(!(xvariable %in% names(para)))
      stop("There is no irregular parameter named", sQuote(xvariable))
    xvalues <- para[[xvariable]]
    xname <- xvariable
  }
  # y variable for plot                  
  if(is.null(coeff)) {
    yvalues <- x$prof
    ylab <- "log PL"
  } else {
    stopifnot(is.character(coeff))
    allcoef <- x$allcoef
    if(!(coeff %in% names(allcoef)))
      stop(paste("There is no coefficient named", sQuote(coeff), "in the fitted model"))
    yvalues <- allcoef[[coeff]]
    ylab <- paste("coefficient:", coeff)
  }
  # start plot
  if(!add)
    do.call.matched("plot.default",
                  resolve.defaults(list(x=range(xvalues), y=range(yvalues)),
                                   list(type="n", main=main),
                                   list(...),
                                   list(ylab=ylab, xlab=xname)))
  # single curve
  if(npara == 1) {
    do.call.matched("lines.default", list(x=xvalues, y=yvalues, ...))
    abline(v = xvalues[x$iopt], lty=2, col="green")
    return(invisible(NULL))
  }

  # multiple curves
  other <- para[, -1, drop=FALSE]
  tapply(1:nrow(para),
         as.list(other),
         function(z, xvalues, yvalues, other, tag) {
           fz <- xvalues[z]
           pz<- yvalues[z]
           lines(fz, pz)
           if(tag) {
             oz <- other[z, , drop=FALSE]
             uniques <- apply(oz, 2, unique)
             labels <- paste(names(uniques), "=", uniques, sep="")
             label <- paste(labels, sep=",")
             ii <- which.max(pz)
             text(fz[ii], pz[ii], label)
           }
         },
         xvalues=xvalues, yvalues=yvalues, other=other, tag=tag)

  abline(v = xvalues[x$iopt], lty=2, col="green")
  return(invisible(NULL))
}
