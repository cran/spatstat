#
#   resolve.defaults.R
#
#  $Revision: 1.3 $ $Date: 2006/06/05 09:29:35 $
#
# Resolve conflicts between several sets of defaults
# Usage:
#     resolve.defaults(list1, list2, list3, .......)
# where the earlier lists have priority 
#
resolve.defaults <- function(...) {
  arglist <- list(...)
  argue <- list()
  if((n <- length(arglist)) > 0)  {
    for(i in seq(n))
      argue <- append(argue, arglist[[i]])
  }
  if(!is.null(nam <- names(argue))) {
    named <- (nam != "")
    arg.unnamed <- argue[!named]
    arg.named <-   argue[named]
    if(any(discard <- duplicated(names(arg.named)))) 
      arg.named <- arg.named[!discard]
    argue <- append(arg.unnamed, arg.named)
  }
  return(argue)
}


do.call.matched <- function(fname, arglist, funargs, extrargs=NULL) {
  fun <- get(fname)
  if(!is.function(fun))
    stop(paste("internal error: function \"", fname, "\" not found",
               sep=""))
  if(missing(funargs))
    funargs <- names(formals(fun))
  funargs <- c(funargs, extrargs)
  givenargs <- names(arglist)
  matched <- givenargs %in% funargs
  do.call(fname, arglist[matched])
}

  
