#
#
#    fv.R
#
#    class "fv" of function value objects
#
#    $Revision: 1.55 $   $Date: 2010/04/16 10:46:57 $
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
  if(!missing(yexp)) {
    if(is.null(yexp)) yexp <- ylab
    else stopifnot(is.language(yexp))
  }
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
    alim <- range(argue[is.finite(argue)], na.rm=TRUE)
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
  attr(x, "dotnames") <- NULL
  # 
  class(x) <- c("fv", class(x))
  return(x)
}

.Spatstat.FvAttrib <- c(
                        "argu",
                        "valu",
                        "ylab",
                        "yexp",
                        "fmla",
                        "alim",
                        "labl",
                        "desc",
                        "units",
                        "fname",
                        "dotnames")

as.data.frame.fv <- function(x, ...) {
  stopifnot(is.fv(x))
  fva <- .Spatstat.FvAttrib
  attributes(x)[fva] <- NULL
  class(x) <- "data.frame"
  x
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

# manipulating the names in fv objects

.Spatstat.FvAbbrev <- c(
                        ".x",
                        ".y",
                        ".",
                        "..")

fvnames <- function(X, a=".") {
  verifyclass(X, "fv")
  if(!is.character(a) || length(a) > 1)
    stop("argument a must be a character string")
  switch(a,
         ".y"={
           return(attr(X, "valu"))
         },
         ".x"={
           return(attr(X, "argu"))
         },
         "." = {
           # The specified 'dotnames'
           dn <- attr(X, "dotnames")
           if(is.null(dn)) 
             dn <- fvnames(X, "..")
           return(dn)
         },
         ".."={
           # all column names other than the function argument
           allvars <- names(X)
           argu <- attr(X, "argu")
           nam <- allvars[allvars != argu]
           nam <- rev(nam) # convention
           return(nam)
         },
         stop(paste("Unrecognised abbreviation", dQuote(a)))
       )
}

"fvnames<-" <- function(X, a=".", value) {
  verifyclass(X, "fv")
  if(!is.character(a) || length(a) > 1)
    stop(paste("argument", sQuote("a"), "must be a character string"))
  if(a == "..") {
    warning(paste("Cannot reset fvnames(x,", dQuote(".."), ")"))
    return(X)
  }
  if(a == "." && length(value) == 0) {
    # clear the dotnames
    attr(X, "dotnames") <- NULL
    return(X)
  }
  # validate the names
  switch(a,
         ".x"=,
         ".y"={
           if(!is.character(value) || length(value) != 1)
             stop("value should be a single string")
         },
         "."={
           if(!is.character(value))
             stop("value should be a character vector")
         },
         stop(paste("Unrecognised abbreviation", dQuote(a)))
       )
  # check the names match existing column names
  tags <- names(X)
  if(any(nbg <- !(value %in% tags))) 
    stop(paste(ngettext(sum(nbg), "The string", "The strings"),
               commasep(dQuote(value[nbg])),
               ngettext(sum(nbg),
                        "does not match the name of any column of X", 
                        "do not match the names of any columns of X")))
  # reassign names
  switch(a,
         ".x"={
           attr(X, "argu") <- value
         },
         ".y"={
           attr(X, "valu") <- value
         },
         "."={
           attr(X, "dotnames") <- value
         })
  return(X)
}

bind.fv <- function(x, y, labl=NULL, desc=NULL, preferred=NULL) {
  verifyclass(x, "fv")
  a <- attributes(x)
  if(is.fv(y)) {
    # y is already an fv object
    # check compatibility of 'r' values
    xr <- attr(x, "argu")
    yr <- attr(y, "argu")
    rx <- x[[xr]]
    ry <- y[[yr]]
    if((length(rx) != length(rx)) || 
               (max(abs(rx-ry)) > .Machine$double.eps))
      stop("fv objects x and y have incompatible domains")
    # reduce y to data frame and strip off 'r' values
    ystrip <- as.data.frame(y)
    yrpos <- which(colnames(ystrip) == yr)
    ystrip <- ystrip[, -yrpos, drop=FALSE]
    # determine descriptors
    if(is.null(labl)) labl <- attr(y, "labl")[-yrpos]
    if(is.null(desc)) desc <- attr(y, "desc")[-yrpos]
    #
    y <- ystrip
  } else {
    # y is a matrix or data frame
    y <- as.data.frame(y)
  }
  
  # check for duplicated column names
  allnames <- c(colnames(x), colnames(y))
  if(any(dup <- duplicated(allnames))) {
    nbg <- unique(allnames[dup])
    nn <- length(nbg)
    warning(paste("The column",
                  ngettext(nn, "name", "names"),
                  commasep(sQuote(nbg)),
                  ngettext(nn, "was", "were"),
                  "duplicated. Unique names were generated"))
    allnames <- make.names(allnames, unique=TRUE, allow_ = FALSE)
    colnames(y) <- allnames[ncol(x) + seq(ncol(y))]
  }
      
  if(is.null(labl))
    labl <- paste("%s[", colnames(y), "](r)", sep="")
  else if(length(labl) != ncol(y))
    stop(paste("length of", sQuote("labl"),
               "does not match number of columns of y"))
  if(is.null(desc))
    desc <- character(ncol(y))
  else if(length(desc) != ncol(y))
    stop(paste("length of", sQuote("desc"),
               "does not match number of columns of y"))
  if(is.null(preferred))
    preferred <- a$valu

  xy <- cbind(as.data.frame(x), y)
  z <- fv(xy, a$argu, a$ylab, preferred, a$fmla, a$alim,
          c(attr(x, "labl"), labl),
          c(attr(x, "desc"), desc),
          unitname=unitname(a),
          fname=attr(x, "fname"))
  return(z)
}

cbind.fv <- function(...) {
  a <- list(...)
  n <- length(a)
  if(n == 0)
    return(NULL)
  if(n == 1) {
    # single argument - extract it
    a <- a[[1]]
    # could be an fv object 
    if(is.fv(a))
      return(a)
    n <- length(a)
  }
  z <- a[[1]]
  if(!is.fv(z))
    stop("First argument should be an object of class fv")
  if(n > 1)
    for(i in 2:n) 
      z <- bind.fv(z, a[[i]])
  return(z)
}

collapse.fv <- function(..., same=NULL, different=NULL) {
  x <- list(...)
  n <- length(x)
  if(n == 0)
    return(NULL)
  if(n == 1)  {
    # single argument - could be a list - extract it
    x1 <- x[[1]]
    if(!is.fv(x1))
      x <- x1
  } 
  if(!all(unlist(lapply(x, is.fv))))
    stop("arguments should be objects of class fv")
  if(is.null(same)) same <- character(0)
  if(is.null(different)) different <- character(0)
  if(any(duplicated(c(same, different))))
    stop(paste("The arguments", sQuote("same"), "and", sQuote("different"),
               "should not have entries in common"))
  either <- c(same, different)
  # names for different versions
  versionnames <- names(x)
  if(is.null(versionnames))
    versionnames <- paste("x", seq(length(x)), sep="")
  shortnames <- abbreviate(versionnames)
  # extract the common values
  y <- x[[1]]
  if(!(fvnames(y, ".y") %in% same))
    fvnames(y, ".y") <- same[1]
  z <- y[, c(fvnames(y, ".x"), same)]
  dotnames <- same
  # now merge the different values
  for(i in seq(length(x))) {
    # extract values for i-th object
    xi <- x[[i]]
    wanted <- (names(xi) %in% different)
    y <- as.data.frame(xi)[, wanted, drop=FALSE]
    desc <- attr(xi, "desc")[wanted]
    labl <- attr(xi, "labl")[wanted]
    # relabel
    prefix <- shortnames[i]
    preamble <- versionnames[i]
    names(y) <- if(ncol(y) == 1) prefix else paste(prefix,names(y),sep="")
    dotnames <- c(dotnames, names(y))
    # glue onto fv object
    z <- bind.fv(z, y,
                 labl=paste(prefix, labl, sep=""),
                 desc=paste(preamble, desc))
  }
  fvnames(z, ".") <- dotnames
  return(z)
}


# rename one of the columns of an fv object
tweak.fv.entry <- function(x, current.tag, new.labl=NULL, new.desc=NULL, new.tag=NULL) {
  hit <- (names(x) == current.tag)
  if(!any(hit))
    return(x)
  # update descriptions of column
  i <- min(which(hit))
  if(!is.null(new.labl)) attr(x, "labl")[i] <- new.labl
  if(!is.null(new.desc)) attr(x, "desc")[i] <- new.desc
  # adjust column tag
  if(!is.null(new.tag)) {
    names(x)[i] <- new.tag
    # update dotnames
    dn <- fvnames(x, ".")
    if(current.tag %in% dn ) {
      dn[dn == current.tag] <- new.tag
      fvnames(x, ".") <- dn
    }
    # if the tweaked column is the preferred value, adjust accordingly
    if(attr(x, "valu") == current.tag)
      attr(x, "valu") <- new.tag
    # if the tweaked column is the function argument, adjust accordingly
    if(attr(x, "argu") == current.tag)
      attr(x, "valu") <- new.tag
  }
  return(x)
}


# change some or all of the auxiliary text in an fv object
rebadge.fv <- function(x, new.ylab, new.fname,
                       tags, new.desc, new.labl,
                       new.yexp=new.ylab, new.dotnames,
                       new.preferred, new.formula) {
  if(!missing(new.ylab)) 
    attr(x, "ylab") <- new.ylab
  if(!missing(new.yexp) || !missing(new.ylab))
    attr(x, "yexp") <- new.yexp
  if(!missing(new.fname))
    attr(x, "fname") <- new.fname
  if(!missing(tags) && !(missing(new.desc) && missing(new.labl))) {
    nama <- names(x)
    desc <- attr(x, "desc")
    labl <- attr(x, "labl")
    for(i in seq(length(tags)))
    if(!is.na(m <- match(tags[i], nama))) {
      if(!missing(new.desc)) desc[m] <- new.desc[i]
      if(!missing(new.labl)) labl[m] <- new.labl[i]
    }
    attr(x, "desc") <- desc
    attr(x, "labl") <- labl
  }
  if(!missing(new.dotnames))
    fvnames(x, ".") <- new.dotnames
  if(!missing(new.preferred)) {
    stopifnot(new.preferred %in% names(x))
    attr(x, "valu") <- new.preferred
  }
  if(!missing(new.formula))
    attr(x, "fmla") <- new.formula
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
  dnames <- fvnames(data, ".")
  xname <- fvnames(data, ".x")
  yname <- fvnames(data, ".y")
  ud <- as.call(lapply(c("cbind", dnames), as.name))
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


  

