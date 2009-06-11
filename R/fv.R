#
#
#    fv.R
#
#    class "fv" of function value objects
#
#    $Revision: 1.39 $   $Date: 2009/06/02 23:18:34 $
#
#
#    An "fv" object represents one or more related functions
#    of the same argument, such as different estimates of the K function.
#
#    It is a data.frame with additional attributes
#    
#         argu       column name of the function argument (typically "r")
#
#         valu       column name of the recommended function
#
#         ylab       generic label for y axis e.g. K(r)
#
#         fmla       default plot formula
#
#         alim       recommended range of function argument
#
#         labl       recommended xlab/ylab for each column
#
#         desc       longer description for each column 
#
#         unitname   name of unit of length for 'r'
#
#    Objects of this class are returned by Kest(), etc
#
##################################################################
# creator

fv <- function(x, argu="r", ylab=NULL, valu, fmla=NULL,
               alim=NULL, labl=names(x), desc=NULL, unitname=NULL,
               fname=NULL, yexp=ylab) {
  stopifnot(is.data.frame(x))
  # check arguments
  stopifnot(is.character(argu))
  if(!is.null(ylab))
    stopifnot(is.character(ylab) || is.language(ylab))
  if(!missing(yexp))
    stopifnot(is.language(yexp))
  stopifnot(is.character(valu))
  
  if(!(argu %in% names(x)))
    stop(paste(sQuote("argu"), "must be the name of a column of x"))

  if(!(valu %in% names(x)))
    stop(paste(sQuote("valu"), "must be the name of a column of x"))

  if(is.null(fmla))
    fmla <- as.formula(paste(valu, "~", argu))
  else if(!inherits(fmla, "formula") && !is.character(fmla))
    stop(paste(sQuote("fmla"), "should be a formula or a string"))
  # convert to string
  fmla <- deparse(fmla)

  if(is.null(alim)) {
    argue <- x[[argu]]
    xlim <- range(argue[is.finite(argue)], na.rm=TRUE)
  }
  if(!is.numeric(alim) || length(alim) != 2)
    stop(paste(sQuote("alim"), "should be a vector of length 2"))
  if(!is.character(labl))
    stop(paste(sQuote("labl"), "should be a vector of strings"))
  stopifnot(length(labl) == ncol(x))
  if(is.null(desc))
    desc <- character(ncol(x))
  else {
    stopifnot(is.character(desc))
    stopifnot(length(desc) == ncol(x))
    nbg <- is.na(desc)
    if(any(nbg)) desc[nbg] <- ""
  }
  if(!is.null(fname))
    stopifnot(is.character(fname) && length(fname) == 1)
  # pack attributes
  attr(x, "argu") <- argu
  attr(x, "valu") <- valu
  attr(x, "ylab") <- ylab
  attr(x, "yexp") <- yexp
  attr(x, "fmla") <- fmla
  attr(x, "alim") <- alim
  attr(x, "labl") <- labl
  attr(x, "desc") <- desc
  attr(x, "units") <- as.units(unitname)
  attr(x, "fname") <- fname
  # 
  class(x) <- c("fv", class(x))
  return(x)
}

is.fv <- function(x) {
  inherits(x, "fv")
}

as.fv <- function(x) {
  if(is.fv(x))
    return(x)
  else if(inherits(x, "data.frame"))
    return(fv(x, names(x)[1], , names(x)[2]))
  else if(inherits(x, "fasp") && length(x$which) == 1)
    return(x$fns[[1]])
  else
    stop(paste("Don't know how to convert this to an object of class",
               sQuote("fv")))
}

print.fv <- function(x, ...) {
  verifyclass(x, "fv")
  nama <- names(x)
  a <- attributes(x)
  cat(paste("Function value object (class ", sQuote("fv"), ")\n", sep=""))
  if(!is.null(ylab <- a$ylab)) {
    if(is.language(ylab))
      ylab <- deparse(ylab)
    cat(paste("for the function", a$argu, "->", ylab, "\n"))
  }
  # Descriptions ..
  desc <- a$desc
  # .. may require insertion of ylab
  if(!is.null(ylab))
    desc <- sprintf(desc, ylab)
  # Labels ..
  labl <- a$labl
  # .. may require insertion of function name if it is known
  if(!is.null(fname <- attr(x, "fname")))
    labl <- sprintf(labl, fname)
  # Start printing
  cat("Entries:\n")
  lablen <- nchar(labl)
  labjump <- max(c(lablen,5)) + 3
  idlen <- nchar(nama)
  idjump <- max(c(idlen,5)) + 3
  pad <- function(n) { paste(rep(" ", n), collapse="") }
  cat("id", pad(idjump-2), "label", pad(labjump - 5), "description\n", sep="")
  cat("--", pad(idjump-2), "-----", pad(labjump - 5), "-----------\n", sep="")
  for(j in seq(ncol(x))) 
    cat(paste(nama[j], pad(idjump - idlen[j]),
              labl[j],pad(labjump - lablen[j]),
              desc[j],"\n", sep=""))
  cat("--------------------------------------\n\n")
  cat("Default plot formula:\n\t")
  print.formula(as.formula(a$fmla))
  alim <- signif(a$alim, 5)
  cat(paste("\nRecommended range of argument ", a$argu,
            ": [", alim[1], ", ", alim[2], "]\n", sep=""))
  ledge <- summary(unitname(x))$legend
  if(!is.null(ledge))
    cat(paste(ledge, "\n"))
  invisible(NULL)
}

bind.fv <- function(x, y, labl, desc, preferred) {
  verifyclass(x, "fv")
  y <- as.data.frame(y)
  a <- attributes(x)
  
  if(length(labl) != ncol(y))
    stop(paste("length of", sQuote("labl"),
               "does not match number of columns of y"))
  if(missing(desc) || is.null(desc))
    desc <- character(ncol(y))
  else if(length(desc) != ncol(y))
    stop(paste("length of", sQuote("desc"),
               "does not match number of columns of y"))
  if(missing(preferred))
    preferred <- a$valu

  xy <- cbind(as.data.frame(x), y)
  z <- fv(xy, a$argu, a$ylab, preferred, a$fmla, a$alim,
          c(attr(x, "labl"), labl),
          c(attr(x, "desc"), desc),
          unitname=unitname(a),
          fname=attr(x, "fname"))
  return(z)
}

# rename one of the columns of an fv object
tweak.fv.entry <- function(x, current.tag, new.labl=NULL, new.desc=NULL) {
  hit <- (names(x) == current.tag)
  if(any(hit)) {
    i <- min(which(hit))
    if(!is.null(new.labl)) attr(x, "labl")[i] <- new.labl
    if(!is.null(new.desc)) attr(x, "desc")[i] <- new.desc
  }
  return(x)
}

# change all the text in an fv object
rebadge.fv <- function(x, new.ylab=NULL, new.fname=NULL,
                       tags=NULL, new.desc=NULL, new.labl=NULL,
                       new.yexp=new.ylab) {
  if(!is.null(new.ylab))
    attr(x, "ylab") <- new.ylab
  if(!is.null(new.yexp))
    attr(x, "yexp") <- new.yexp
  if(!missing(new.fname))
    attr(x, "fname") <- new.fname
  if(!is.null(tags) && !(is.null(new.desc) && is.null(new.labl))) {
    nama <- names(x)
    desc <- attr(x, "desc")
    labl <- attr(x, "labl")
    for(i in seq(length(tags)))
    if(!is.na(m <- match(tags[i], nama))) {
      if(!is.null(new.desc)) desc[m] <- new.desc[i]
      if(!is.null(new.labl)) labl[m] <- new.labl[i]
    }
    attr(x, "desc") <- desc
    attr(x, "labl") <- labl
  }
  return(x)
}

# subset extraction operator
"[.fv" <-
  function(x, i, j, ..., drop=FALSE)
{
  Nindices <- (!missing(i)) + (!missing(j))
  if(Nindices == 0)
    return(x)
  y <- as.data.frame(x)
  if(Nindices == 2)
    z <- y[i, j, drop=FALSE]
  else if(!missing(i))
    z <- y[i, , drop=FALSE]
  else
    z <- y[ , j, drop=FALSE]

  if(missing(j)) 
    selected <- seq(ncol(x))
  else {
    nameindices <- seq(names(x))
    names(nameindices) <- names(x)
    selected <- as.vector(nameindices[j])
  }

  nama <- names(z)
  argu <- attr(x, "argu")
  if(!(argu %in% nama))
    stop(paste("The function argument", sQuote(argu), "must not be removed"))
  valu <- attr(x, "valu")
  if(!(valu %in% nama))
    stop(paste("The default column of function values",
               sQuote(valu), "must not be removed"))

  # If range of argument was implicitly changed, adjust "alim"
  alim <- attr(x, "alim")
  rang <- range(z[[argu]])
  alim <- c(max(alim[1], rang[1]),
            min(alim[2], rang[2]))

  return(fv(z, argu=attr(x, "argu"),
               ylab=attr(x, "ylab"),
               valu=attr(x, "valu"),
               fmla=attr(x, "fmla"),
               alim=alim,
               labl=attr(x, "labl")[selected],
               desc=attr(x, "desc")[selected],
               unitname=attr(x, "units"),
               fname=attr(x,"fname")))
}  

#   method for with()

with.fv <- function(data, expr, ..., drop=TRUE) {
  verifyclass(data, "fv")
  # convert syntactic expression to 'expression' object
  e <- as.expression(substitute(expr))
  # convert syntactic expression to call
  elang <- substitute(expr)
  # expand "."
  dotnames <- fvnames(data, ".")
  xname <- fvnames(data, ".x")
  yname <- fvnames(data, ".y")
  ud <- as.call(lapply(c("cbind", dotnames), as.name))
  ux <- as.name(xname)
  uy <- as.name(yname)
  elang <- eval(substitute(substitute(ee,
                                      list(.=ud, .x=ux, .y=uy)),
                           list(ee=elang)))
  # evaluate expression
  datadf <- as.data.frame(data)
  results <- eval(elang, as.list(datadf))
  # --------------------
  # make sense of the results
  #
  nx <- nrow(datadf)
  # 
  if(!is.matrix(results) && !is.data.frame(results)) {
    # result is a vector
    if(length(results) != nx) {
      # format not understood
#      warning("Calculation produced a vector of the wrong length")
      return(results)
    }
    # result is a vector of the right length
    if(drop)
      return(as.vector(results))
    else
      results <- matrix(results, nrow=nx, ncol=1)
  }
  # result is a matrix or data frame
  if(nrow(results) != nx) {
    # format not understood - dump the values
#    warning("Calculation yielded a matrix or data frame of the wrong dimensions")
    return(results)
  }
  # result is a matrix or data frame of the right dimensions
  # make a new fv object
  # ensure columns of results have names
  if(is.null(colnames(results)))
    colnames(results) <- paste("col", seq(ncol(results)), sep="")
  resultnames <- colnames(results)
  # get values of function argument
  xvalues <- datadf[[xname]]
  # tack onto result matrix
  results <- cbind(xvalues, results)
  colnames(results) <- c(xname, resultnames)
  results <- data.frame(results)
  # check for alteration of column names
  oldnames <- resultnames
  resultnames <- colnames(results)[-1]
  if(any(resultnames != oldnames))
    warning("some column names were illegal and have been changed")
  # Build up fv object
  # decide which of the columns should be the preferred value
  newyname <- if(yname %in% resultnames) yname else resultnames[1]
  # construct default plot formula
  lhs <- resultnames
  if(length(lhs) > 1)
    lhs <- paste("cbind", paren(paste(lhs, collapse=", ")))
  fmla <- as.formula(paste(lhs, "~", xname))
  # construct description strings
  desc <- c(attr(data, "desc")[1], paste("Computed value", resultnames))
  # form fv object and return
  out <- fv(results, argu=xname, valu=newyname,
            desc=desc, alim=attr(data, "alim"), fmla=fmla, 
            unitname=unitname(data), fname="?")
  return(out)
}


# translate obscure names

fvnames <- function(X, a=".") {
  verifyclass(X, "fv")
  if(!is.character(a) || length(a) > 1)
    stop("argument a must be a character string")
  switch(a,
         "." = {
           dn <- attr(X, "dotnames")
           if(is.null(dn)) {
             argu <- attr(X, "argu")
             allvars <- names(X)
             dn <- allvars[allvars != argu]
             dn <- rev(dn) # convention
           }
           return(dn)
         },
         ".y"={
           return(attr(X, "valu"))
         },
         ".x"={
           return(attr(X, "argu"))
         },
         stop(paste("Unrecognised abbreviation", dQuote(a)))
       )
}

  
# stieltjes integration for fv objects

stieltjes <- function(f, M, ...) {
  # stieltjes integral of f(x) dM(x)
  if(!is.fv(M))
    stop("M must be an object of class fv")
  if(!is.function(f))
    stop("f must be a function")
  # integration variable
  argu <- attr(M, "argu")
  x <- M[[argu]]
  # values of integrand
  fx <- f(x, ...)
  # estimates of measure
  valuenames <- names(M) [names(M) != argu]
  Mother <- as.data.frame(M)[, valuenames]
  # increments of measure
  dM <- apply(Mother, 2, diff)
  dM <- rbind(dM, 0)
  dM[is.na(dM)] <- 0
  # integrate f(x) dM(x)
  results <- apply(fx * dM, 2, sum)
  return(as.list(results))
}


  

