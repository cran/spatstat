#
# layered.R
#
# Simple mechanism for layered plotting
#
#  $Revision: 1.1 $  $Date: 2011/09/01 04:11:32 $
#

layered <- function(..., plotargs=NULL) {
  out <- list(...)
  if(!is.null(plotargs)) {
    if(!is.list(plotargs) || !all(unlist(lapply(plotargs, is.list))))
      stop("plotargs should be a list of lists")
    if(length(plotargs) != length(out))
      stop("plotargs should have one component for each element of the list")
    names(plotargs) <- names(out)
  }
  attr(out, "plotargs") <- plotargs
  class(out) <- c("layered", class(out))
  return(out)
}

print.layered <- function(x, ...) {
  cat("Layered object\n")
  for(i in seq_along(x)) {
    cat(paste("\nLayer ", i, ":\n", sep=""))
    print(x[[i]])
  }
  invisible(NULL)
}

plot.layered <- function(x, ..., which=NULL, plotargs=NULL) {
  xname <- deparse(substitute(x))
  xp <- if(is.null(which)) x else x[which]
  if(length(xp) == 0)
    return(invisible(NULL))
  if(is.null(plotargs)) {
    plotargs <- attr(x, "plotargs")
    if(!is.null(plotargs) && !is.null(which)) plotargs <- plotargs[which]
  } else {
    if(!is.list(plotargs) || !all(unlist(lapply(plotargs, is.list))))
      stop("plotargs should be a list of lists")
    if(length(plotargs) != length(xp))
      stop("plotargs should have one component for each layer to be plotted")
  }
  out <- list()
  for(i in seq_along(xp)) {
    first <- (i == 1)
    out[[i]] <-
      do.call("plot",
              resolve.defaults(list(x=xp[[i]]),
                               plotargs[[i]],
                               list(...),
                               if(first) list(main=xname) else list(add=TRUE)))
  }
  return(invisible(out))
}

