#
#   resolve.defaults.R
#
#  $Revision: 1.9 $ $Date: 2011/09/06 03:04:34 $
#
# Resolve conflicts between several sets of defaults
# Usage:
#     resolve.defaults(list1, list2, list3, .......)
# where the earlier lists have priority 
#
resolve.defaults <- function(..., .StripNull=FALSE) {
  # Each argument is a list. Append them.
  argue <- c(...)
  if(!is.null(nam <- names(argue))) {
    named <- nzchar(nam)
    arg.unnamed <- argue[!named]
    arg.named <-   argue[named]
    if(any(discard <- duplicated(names(arg.named)))) 
      arg.named <- arg.named[!discard]
    argue <- append(arg.unnamed, arg.named)
  }
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


  
