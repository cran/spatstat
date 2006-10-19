#
#	fasp.R
#
#	$Revision: 1.16 $	$Date: 2006/10/09 03:17:44 $
#
#
#-----------------------------------------------------------------------------
#

# creator
fasp <- function(fns, titles, formulae, which,
                 dataname=NULL, title=NULL) {
  stopifnot(is.matrix(which))
  stopifnot(is.list(fns))
  stopifnot(length(fns) == length(which))

  fns <- lapply(fns, as.fv)
  n   <- length(which)

  if(!missing(titles)) 
    stopifnot(length(titles) == length(which))
  else
    titles <- lapply(seq(length(which)), function(i) NULL)

  if(!missing(formulae)) {
    # verify format and convert to character vector
    formulae <- FormatFaspFormulae(formulae, "formulae")
    # ensure length matches length of "fns"
    if(length(formulae) == 1 && n > 1)
        # single formula - replicate it
        formulae <- rep(formulae, n)
    else 
        stopifnot(length(formulae) == length(which))
  }

  rslt <- list(fns=fns, titles=titles, default.formula=formulae,
               which=which, dataname=dataname, title=title)
  class(rslt) <- "fasp"
  return(rslt)
}

# subset operator

"[.fasp" <-
"subset.fasp" <-
  function(x, I, J, drop, ...) {

        verifyclass(x, "fasp")
        
        m <- nrow(x$which)
        n <- ncol(x$which)
        
        if(missing(I)) I <- 1:m
        if(missing(J)) J <- 1:n

        # determine index subset for lists 'fns', 'titles' etc
        included <- rep(FALSE, length(x$fns))
        w <- as.vector(x$which[I,J])
        if(!any(w))
          stop("result is empty")
        included[w] <- TRUE

        # determine positions in shortened lists
        whichIJ <- x$which[I,J,drop=FALSE]
        newk <- cumsum(included)
        newwhich <- matrix(newk[whichIJ],
                           ncol=ncol(whichIJ), nrow=nrow(whichIJ))

        # create new fasp object
        Y <- fasp(fns      = x$fns[included],
                  titles   = x$titles[included],
                  formulae = x$default.formula[included],
                  which    = newwhich,
                  dataname = x$dataname,
                  title    = x$title)
        return(Y)
}


print.fasp <- function(x, ...) {
  verifyclass(x, "fasp")
  cat(paste("Function array (class", sQuote("fasp"), "\n"))
  dim <- dim(x$which)
  cat(paste("Dimensions: ", dim[1], "x", dim[2], "\n"))
  cat(paste("Title:", if(is.null(x$title)) "(None)" else x$title, "\n"))
  invisible(NULL)
}

FormatFaspFormulae <- function(f, argname) {
  # f should be a single formula object, a list of formula objects,
  # a character vector, or a list containing formulae and strings.
  # It will be converted to a character vector.
  
  zapit <- function(x, argname) {
    if(inherits(x, "formula")) deparse(x)
    else if(is.character(x)) x
    else stop(paste("The entries of",
                    sQuote(argname),
                    "must be formula objects or strings"))
  }

  result <-
    if(is.character(f))
      f
    else if(inherits(f, "formula"))
      deparse(f)
    else if(is.list(f))
      unlist(lapply(f, zapit, argname=argname))
    else stop(paste(sQuote(argname),
                    "should be a formula, a list of formulae,",
                    "or a character vector"))

  return(result)
}
