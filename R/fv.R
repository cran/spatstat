#
#
#    fv.R
#
#    class "fv" of function value objects
#
#    $Revision: 1.14 $   $Date: 2006/04/12 18:00:04 $
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
#    Objects of this class are returned by Kest(), etc
#
##################################################################
# creator

fv <- function(x, argu="r", ylab=NULL, valu, fmla=NULL,
               alim=NULL, labl=names(x), desc=NULL) {
  stopifnot(is.data.frame(x))
  # check arguments
  stopifnot(is.character(argu))
  if(!is.null(ylab))
    stopifnot(is.character(ylab) || is.language(ylab))
  stopifnot(is.character(valu))
  
  if(!(argu %in% names(x)))
    stop("\`argu\' must be the name of a column of x")

  if(!(valu %in% names(x)))
    stop("\`valu\' must be the name of a column of x")

  if(is.null(fmla))
    fmla <- as.formula(paste(valu, "~", argu))
  else if(!inherits(fmla, "formula") && !is.character(fmla))
    stop("\`fmla\' should be a formula or a string")
  # convert to string
  fmla <- deparse(fmla)

  if(is.null(alim)) {
    argue <- x[[argu]]
    xlim <- range(argue[is.finite(argue)], na.rm=TRUE)
  }
  if(!is.numeric(alim) || length(alim) != 2)
    stop("\`alim\' should be a vector of length 2")
  if(!is.character(labl))
    stop("\`labl\' should be a vector of strings")
  stopifnot(length(labl) == ncol(x))
  if(is.null(desc))
    desc <- character(ncol(x))
  else {
    stopifnot(is.character(desc))
    stopifnot(length(desc) == ncol(x))
    nbg <- is.na(desc)
    if(any(nbg)) desc[nbg] <- ""
  }
  # pack attributes
  attr(x, "argu") <- argu
  attr(x, "valu") <- valu
  attr(x, "ylab") <- ylab
  attr(x, "fmla") <- fmla
  attr(x, "alim") <- alim
  attr(x, "labl") <- labl
  attr(x, "desc") <- desc
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
  else if(inherits(x, "fasp") && length(which) == 1)
    return(x$funs[[1]])
  else
    stop("Don't know how to convert this to an \"fv\" object")
}

print.fv <- function(x, ...) {
  verifyclass(x, "fv")
  nama <- names(x)
  a <- attributes(x)
  cat("Function value object (class \"fv\")\n")
  if(!is.null(ylab <- a$ylab)) {
    if(is.language(ylab))
      ylab <- deparse(ylab)
    cat(paste("for the function", a$argu, "->", ylab, "\n"))
  }
  cat("Entries:\n")
  len <- nchar(a$labl)
  tabjump <- max(c(len, 5)) + 3
  pad <- function(n) { paste(character(n),collapse=" ") }
  cat("id\tlabel", pad(tabjump - 5), "description\n", sep="")
  cat("--\t-----", pad(tabjump - 5), "-----------\n", sep="")
  for(j in seq(ncol(x))) 
    cat(paste(nama[j],"\t",
              a$labl[j],pad(tabjump - len[j]),
              a$desc[j],"\n", sep=""))
  cat("--------------------------------------\n\n")
  cat("Default plot formula:\n\t")
  print.formula(as.formula(a$fmla))
  cat(paste("\nRecommended range of argument ", a$argu,
            ": [", a$alim[1], ", ", a$alim[2], "]\n", sep=""))
  invisible(NULL)
}

bind.fv <- function(x, y, labl, desc, preferred) {
  verifyclass(x, "fv")
  y <- as.data.frame(y)
  a <- attributes(x)
  
  if(length(labl) != ncol(y))
    stop("length of \`labl\' does not match number of columns of y")
  if(missing(desc) || is.null(desc))
    desc <- character(ncol(y))
  else if(length(desc) != ncol(y))
    stop("length of \`desc\' does not match number of columns of y")
  if(missing(preferred))
    preferred <- a$valu

  xy <- cbind(as.data.frame(x), y)
  z <- fv(xy, a$argu, a$ylab, preferred, a$fmla, a$alim,
          c(attr(x, "labl"), labl),
          c(attr(x, "desc"), desc))
  return(z)
}

"[.fv" <- subset.fv <- function(x, i, j, ..., drop=FALSE)
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
    stop(paste("The function argument \`", argu, "\' must not be removed",
               sep=""))
  valu <- attr(x, "valu")
  if(!(valu %in% nama))
    stop(paste("The default column of function values \'", valu,
                  "\' must not be removed"))

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
               desc=attr(x, "desc")[selected]))
}  


