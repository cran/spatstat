#
#	fasp.R
#
#	$Revision: 1.10 $	$Date: 2004/01/13 09:56:35 $
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

  if(!missing(titles)) 
    stopifnot(length(titles) == length(which))
  else
    titles <- lapply(seq(length(which)), function(i) NULL)

  if(!missing(formulae)) {
    if(inherits(formulae, "formula")) 
      # single formula 
      # make it a list of same length as "fns"
      formulae <- lapply(seq(length(which)),
                                 function(i, f) {f}, f=formulae)
    else {
      # list of formulae
      stopifnot(is.list(formulae))
      if(!all(unlist(lapply(formulae, inherits, what="formula"))))
        stop("The entries of \`formulae\' should all be formulae")
      if(length(formulae) == 1 && length(which) != 1)
        # list of length 1 - replicate
        formulae <- lapply(seq(length(which)),
                                 function(i, f) {f}, f=formulae[[1]])
      else 
        stopifnot(length(formulae) == length(which))
    }
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
  cat("Function array (class \"fasp\")\n")
  dim <- dim(x$which)
  cat(paste("Dimensions: ", dim[1], "x", dim[2], "\n"))
  cat(paste("Title:", if(is.null(x$title)) "(None)" else x$title, "\n"))
  invisible(NULL)
}
