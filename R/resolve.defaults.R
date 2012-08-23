#
#   resolve.defaults.R
#
#  $Revision: 1.10 $ $Date: 2012/08/21 05:21:09 $
#
# Resolve conflicts between several sets of defaults
# Usage:
#     resolve.defaults(list1, list2, list3, .......)
# where the earlier lists have priority 
#
resolve.defaults <- function(..., .MatchNull=TRUE, .StripNull=FALSE) {
  # Each argument is a list. Append them.
  argue <- c(...)
  # is NULL a possible value?
  if(!.MatchNull) {
    isnul <- unlist(lapply(argue, is.null))
    argue <- argue[!isnul]
  }
  if(!is.null(nam <- names(argue))) {
    named <- nzchar(nam)
    arg.unnamed <- argue[!named]
    arg.named <-   argue[named]
    if(any(discard <- duplicated(names(arg.named)))) 
      arg.named <- arg.named[!discard]
    argue <- append(arg.unnamed, arg.named)
  }
  # should NULL become a missing argument?
  if(.StripNull) {
    isnull <- sapply(argue, is.null)
    argue <- argue[!isnull]
  }
  return(argue)
}

do.call.matched <- function(fun, arglist, funargs, extrargs=NULL) {
  if(!is.function(fun) && !is.character(fun))
    stop("Internal error: wrong argument type in do.call.matched")
  if(is.character(fun)) {
    fname <- fun
    fun <- get(fname, mode="function")
    if(!is.function(fun))
      stop(paste("internal error: function", sQuote(fname), "not found",
                 sep=""))
  } 
  if(missing(funargs))
    funargs <- names(formals(fun))
  funargs <- c(funargs, extrargs)
  givenargs <- names(arglist)
  matched <- givenargs %in% funargs
  do.call(fun, arglist[matched])
}

resolve.1.default <- function(.A, ...) {
  res <- resolve.defaults(...)
  hit <- (names(res) == .A)
  if(!any(hit)) return(NULL)
  return(res[[min(which(hit))]])
}
