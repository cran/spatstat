#
# marks.R
#
# stuff for handling marks
#
#

.Spatstat.Forbids.Df <- TRUE

marks <- function(x, ...) {
  UseMethod("marks")
}

# The 'dfok' switch is temporary
# while we convert the code to accept data frames of marks

marks.ppp <- function(x, ..., dfok=FALSE) {
  ma <- x$marks
  if(is.data.frame(ma) && !dfok)
    stop("Sorry, not implemented when the marks are a data frame")
  return(ma)
}

# ------------------------------------------------------------------

"marks<-" <- function(x, ..., value) {
  UseMethod("marks<-")
}

"marks<-.ppp" <- function(x, dfok=FALSE, ..., value) {
  if(is.data.frame(value) && !dfok)
    stop("Sorry, data frames of marks are not yet implemented")
  x <- setmarks(x, value)
  return(x)
}


# -------------------------------------------------

markformat <- function(x) {
  UseMethod("markformat")
}

markformat.ppp <- function(x) {
  mf <- x$markformat
  if(is.null(mf)) {
    x <- as.ppp(x)
    mf <- x$markformat
    if(is.null(mf))
      stop("Internal error: markformat.ppp failed")
  }
  return(mf)
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
  if(any(is.na(marx)))
    switch(na.action,
           warn = {
             warning(paste("some mark values are NA in the point pattern",
                           deparse(substitute(X))), call.=FALSE)
           },
           fatal = {
             return(FALSE)
           },
           ignore = {}
           )
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
  marx <- marks(X, ...)
  if(is.null(marx))
    return(FALSE)
  if(any(is.na(marx)))
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
  return(!is.data.frame(marx) && is.factor(marx))
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
  format <- findmarktype(x)
  switch(format,
         none={return(NULL)},
         vector={ return(c(x,y)) },
         dataframe={ return(rbind(x,y)) },
         listof = {
           z <- append(x,y)
           if(!inherits(z, "listof"))
             z <- as.listof(z)
           return(z)
         },
         stop("Internal error: unrecognised format of marks"))
}

