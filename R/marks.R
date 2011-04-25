#
# marks.R
#
#   $Revision: 1.26 $   $Date: 2011/03/25 04:34:06 $
#
# stuff for handling marks
#
#

marks <- function(x, ...) {
  UseMethod("marks")
}

# The 'dfok' switch is temporary
# while we convert the code to accept data frames of marks

marks.ppp <- function(x, ..., dfok=TRUE) {
  ma <- x$marks
  if((is.data.frame(ma) || is.matrix(ma)) && !dfok)
    stop("Sorry, not implemented when the marks are a data frame")
  return(ma)
}

# ------------------------------------------------------------------

"marks<-" <- function(x, ..., value) {
  UseMethod("marks<-")
}

"marks<-.ppp" <- function(x, ..., dfok=TRUE, value) {
  if((is.data.frame(value) || is.matrix(value)) && !dfok)
    stop("Sorry, data frames of marks are not yet implemented")
  if(is.null(value))
    return(unmark(x))
  m <- value
  if(is.hyperframe(m)) 
    stop("Hyperframes of marks are not supported in ppp objects")
  np <- npoints(x)
  if(!is.data.frame(m) && !is.matrix(m)) {
    # vector of marks
    if(length(m) == 1) m <- rep(m, np)
    else if(np == 0) m <- rep(m, 0) # ensures marked pattern is obtained
    else if(length(m) != np) stop("number of points != number of marks")
    marx <- m
  } else {
    m <- as.data.frame(m)
    # data frame of marks
    if(ncol(m) == 0) {
      # no mark variables
      marx <- NULL
    } else {
      # marks to be attached
      if(nrow(m) == np) {
        marx <- m
      } else {
        # lengths do not match
        if(nrow(m) == 1 || np == 0) {
          # replicate data frame
          marx <- as.data.frame(lapply(as.list(m),
                                       function(x, k) { rep(x, k) },
                                       k=np))
        } else
        stop("number of rows of data frame != number of points")
      }
    }
  }
  # attach/overwrite marks
  Y <- ppp(x$x,x$y,window=x$window,marks=marx, check=FALSE)
  return(Y)
}

"%mark%" <- setmarks <- function(x,value) {
  marks(x) <- value
  return(x)
}

# -------------------------------------------------

markformat <- function(x) {
  UseMethod("markformat")
}

markformat.ppp <- function(x) {
  mf <- x$markformat
  if(is.null(mf)) 
    mf <- markformat(marks(x))
  return(mf)
}

markformat.default <- function(x) {
  if(is.null(x)) return("none")
  if(is.vector(x) || is.factor(x)) return("vector")
  if(is.data.frame(x)) return("dataframe")
  stop("Mark format not understood")
}

# ------------------------------------------------------------------

"is.marked" <-
function(X, ...) {
  UseMethod("is.marked")
}

"is.marked.ppp" <-
function(X, na.action="warn", ...) {
  marx <- marks(X, ...)
  if(is.null(marx))
    return(FALSE)
  if((length(marx) > 0) && any(is.na(marx))) {
    gripe <- paste("some mark values are NA in the point pattern",
                   deparse(substitute(X)))
    switch(na.action,
           warn = warning(gripe, call.=FALSE),
           fatal = stop(gripe, call.=FALSE),
           ignore = {}
           )
  }
  return(TRUE)
}

"is.marked.default" <-
  function(...) { return(FALSE) }


# ------------------------------------------------------------------

is.multitype <- function(X, ...) {
  UseMethod("is.multitype")
}

is.multitype.default <- function(...) { return(FALSE) }

is.multitype.ppp <- function(X, na.action="warn", ...) {
  marx <- marks(X, dfok=TRUE)
  if(is.null(marx))
    return(FALSE)
  if(is.data.frame(marx) && ncol(marx) > 1)
    return(FALSE)
  if(!is.factor(marx))
    return(FALSE)
  if((length(marx) > 0) && any(is.na(marx)))
    switch(na.action,
           warn = {
             warning(paste("some mark values are NA in the point pattern",
                           deparse(substitute(X))))
           },
           fatal = {
             return(FALSE)
           },
           ignore = {}
           )
  return(TRUE)
}

# ------------------------------------------------------------------

unmark <- function(X) {
  UseMethod("unmark")
}

unmark.ppp <- function(X) {
  X$marks <- NULL
  X$markformat <- "none"
  return(X)
}

unmark.splitppp <- function(X) {
  Y <- lapply(X, unmark.ppp)
  class(Y) <- c("splitppp", class(Y))
  return(Y)
}

##### utility functions for subsetting & combining marks #########

findmarktype <- function(x) {
  if(is.null(x)) return("none")
  if(is.vector(x) || is.factor(x)) return("vector")
  if(is.data.frame(x)) return("dataframe")
  if(inherits(x, "listof")) return("listof")
  stop("Internal error: unrecognised mark format")
}

marksubset <- function(x, index, format=NULL) {
  if(is.null(format)) format <- findmarktype(x)
  switch(format,
         none={return(NULL)},
         vector={return(x[index])},
         dataframe={return(x[index,,drop=FALSE])},
         listof={return(x[index])},
         stop("Internal error: unrecognised format of marks"))
}

"%msub%" <- marksubsetop <- function(x,i) { marksubset(x, i) }

"%mrep%" <- markreplicateop <- function(x,n) { 
  format <- findmarktype(x)
  switch(format,
         none={return(NULL)},
         listof=,
         vector={ return(rep(x,n))},
         dataframe={
           return(as.data.frame(lapply(x, function(z, k) rep(z, k), k=n)))
         },
         stop("Internal error: unrecognised format of marks"))
}

"%mapp%" <- markappendop <- function(x,y) { 
  fx <- findmarktype(x)
  fy <- findmarktype(y)
  agree <- (fx == fy)
  if(fx == "data.frame")
    agree <- agree && identical(names(x),names(y)) 
  if(!agree)
    stop("Attempted to concatenate marks that are not compatible")
  switch(fx,
         none   = { return(NULL) },
         vector = {
           if(is.factor(x) || is.factor(y))
             return(cat.factor(x,y))
           else return(c(x,y))
         },
         dataframe = { return(rbind(x,y)) },
         listof = {
           z <- append(x,y)
           if(!inherits(z, "listof"))
             z <- as.listof(z)
           return(z)
         },
         stop("Internal error: unrecognised format of marks"))
}

