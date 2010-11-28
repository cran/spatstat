#
#  rhohat.R
#
#  $Revision: 1.23 $  $Date: 2010/11/24 06:29:12 $
#
#  Non-parametric estimation of a transformation rho(z) determining
#  the intensity function lambda(u) of a point process in terms of a
#  spatial covariate Z(u) through lambda(u) = rho(Z(u)).
#  More generally allows offsets etc.

rhohat <- function(object, covariate, ...,
                   transform=FALSE,
                   n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
                   bwref=bw, covname) {
  if(missing(covname)) 
    covname <- sensiblevarname(deparse(substitute(covariate)), "X")
  callstring <- paste(deparse(sys.call()), collapse = "")  
  # validate model
  if(is.ppp(object) || inherits(object, "quad")) {
    model <- ppm(object, ~1, forcefit=TRUE)
    reference <- "area"
  } else if(is.ppm(object)) {
    model <- object
    reference <- "model"
    modelcall <- model$call
    if(is.null(getglmfit(model)))
      model <- update(model, forcefit=TRUE)
  } else stop("object should be a point pattern or a point process model")

  if(is.null(adjust)) adjust <- 1
  
  # evaluate the covariate at data points and at pixels
  if(is.character(covariate) && length(covariate) == 1) {
    covname <- covariate
    switch(covname,
           x={
             covariate <- function(x,y) { x }
           }, 
           y={
             covariate <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
  }
  stuff <- evalCovar(model, covariate)
  # unpack
  info   <- stuff$info
  values <- stuff$values
  # values at each data point
  ZX      <- values$ZX
  # values at each pixel
  Zvalues <- values$Zvalues
  lambda  <- values$lambda
  # normalising constants
  nX   <- length(ZX)
  area <- area.owin(as.owin(model))
  baseline <- if(reference == "area") area else (mean(lambda) * area)
  kappahat <- nX/baseline
  # estimate densities
  if(is.null(from))
    from <- min(ZX)
  if(is.null(to))
    to   <- max(ZX)
  interpolate <- function(x,y) {
    if(inherits(x, "density") && missing(y))
      approxfun(x$x, x$y, rule=2)
    else 
      approxfun(x, y, rule=2)
  }
  # reference density
  ghat <- density(Zvalues,weights=lambda/sum(lambda),
                  bw=bwref,adjust=adjust,n=n,from=from,to=to, ...)
  xxx <- ghat$x
  ghatfun <- interpolate(ghat)
  # relative density
  if(!transform) {
    # compute ratio of smoothed densities
    fhat <- density(ZX,bw=bw,adjust=adjust,n=n,from=from, to=to, ...)
    fhatfun <- interpolate(fhat)
    yyy <- kappahat * fhatfun(xxx)/ghatfun(xxx)
    # compute variance approximation
    sigma <- fhat$bw
    fstar <- density(ZX,bw=bw,adjust=adjust/sqrt(2),n=n,from=from, to=to, ...)
    fstarfun <- interpolate(fstar)
    const <- 1/(2 * sigma * sqrt(pi))
    vvv  <- const * nX * fstarfun(xxx)/(baseline * ghatfun(xxx))^2
  } else {
    # probability integral transform
    Gfun <- interpolate(ghat$x, cumsum(ghat$y)/sum(ghat$y))
    GZX <- Gfun(ZX)
    # smooth density on [0,1]
    qhat <- density(GZX,bw=bw,adjust=adjust,n=n,from=0, to=1, ...)
    qhatfun <- interpolate(qhat)
    # edge effect correction
    one <- density(seq(0,1,length=512),
                    bw=qhat$bw, adjust=1, n=n,from=0, to=1, ...)
    onefun <- interpolate(one)
    # apply to transformed values
    Gxxx <- Gfun(xxx)
    yyy <- kappahat * qhatfun(Gxxx)/onefun(Gxxx)
    # compute variance approximation
    sigma <- qhat$bw
    qstar <- density(GZX,bw=bw,adjust=adjust/sqrt(2),n=n,from=0, to=1, ...)
    qstarfun <- interpolate(qstar)
    const <- 1/(2 * sigma * sqrt(pi))
    vvv  <- const * nX * qstarfun(Gxxx)/(baseline * onefun(Gxxx))^2
  }
  sd <- sqrt(vvv)
  # pack into fv object
  df <- data.frame(xxx=xxx, rho=yyy, var=vvv, hi=yyy+2*sd, lo=yyy-2*sd)
  names(df)[1] <- covname
  desc <- c(paste("covariate", covname),
            "Estimated covariate effect",
            "Variance of estimator",
            "Upper limit of pointwise 95%% confidence interval",
            "Lower limit of pointwise 95%% confidence interval")
  rslt <- fv(df,
             argu=covname,
             ylab=substitute(rho(X), list(X=as.name(covname))),
             valu="rho",
             fmla= as.formula(paste(". ~ ", covname)),
             alim=c(from,to),
             labl=c(covname,
               paste("%s", paren(covname), sep=""),
               paste("var", paren(covname), sep=""),
               paste("hi", paren(covname), sep=""),
               paste("lo", paren(covname), sep="")),
             desc=desc,
#             unitname=unitname(data.ppm(model)),
             fname="rho",
             yexp=as.expression(substitute(rho(X), list(X=covname))))
  attr(rslt, "dotnames") <- c("rho", "hi", "lo")
  # pack up
  class(rslt) <- c("rhohat", class(rslt))
  attr(rslt,"stuff") <-
    list(reference  = reference,
         modelcall  = if(reference == "model") modelcall else NULL,
         callstring = callstring,
         sigma      = sigma,
         covname    = paste(covname, collapse=""),
         ZX         = ZX,
         transform  = transform)
  return(rslt)
}

print.rhohat <- function(x, ...) {
  s <- attr(x, "stuff")
  cat("Intensity function estimate (class rhohat)\n")
  cat(paste("for the covariate", s$covname, "\n"))
  switch(s$reference,
         area=cat("Function values are absolute intensities\n"),
         model={
           cat("Function values are relative to fitted model\n")
           print(s$modelcall)
         })
  cat("Estimation method: ")
  if(s$transform)
    cat(paste("probability integral transform,",
              "edge-corrected fixed bandwidth kernel smoothing on [0,1]\n"))
  else 
    cat("fixed-bandwidth kernel smoothing\n")
  cat(paste("Call:", s$callstring, "\n"))
  cat(paste("Actual smoothing bandwidth sigma = ", signif(s$sigma,5), "\n"))
  NextMethod("print")
}

plot.rhohat <- function(x, ..., do.rug=TRUE) {
  xname <- deparse(substitute(x))
  s <- attr(x, "stuff")
  covname <- x$covname
  asked.rug <- !missing(do.rug) && identical(rug, TRUE)
  do.call("plot.fv", resolve.defaults(list(x), list(...),
                                      list(main=xname, shade=c("hi", "lo"))))
  if(do.rug) {
    rugx <- ZX <- s$ZX
    # check whether it's the default plot
    argh <- list(...)
    isfo <- unlist(lapply(argh, inherits, what="formula"))
    if(any(isfo)) {
      fmla <- argh[[min(which(isfo))]]
      rhs <- rhs.of.formula(fmla)
      vars <- variablesinformula(rhs)
      vars <- vars[vars %in% names(x)]
      if(length(vars) == 1 && vars == covname) {
        # expression in terms of covariate
        rhstr <- as.character(rhs)[2]
        dat <- list(ZX)
        names(dat) <- covname
        rugx <- as.numeric(eval(parse(text=rhstr), dat))
      } else {
        warning("Unable to add rug plot")
        rugx <- NULL
      }
    } 
    if(!is.null(rugx))
      rug(rugx)
  }
  invisible(NULL)
}

