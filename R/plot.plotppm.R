#
# plot.plotppm.R
#
# engine of plot method for ppm
#
# $Revision: 1.5 $  $Date: 2005/06/12 20:31:20 $
#
#

plot.plotppm <- function(x,data=NULL,trend=TRUE,cif=TRUE,pause=TRUE,
                         how=c("persp","image","contour"), ...)
{
  verifyclass(x,"plotppm")
  
  # determine main plotting actions
  superimpose <- !is.null(data)
  if(!missing(trend) && (trend & is.null(x[["trend"]])))
    stop("No trend to plot.\n")
  trend <- trend & !is.null(x[["trend"]])
  if(!missing(cif) && (cif & is.null(x[["cif"]])))
    stop("No cif to plot.\n")
  cif <- cif & !is.null(x[["cif"]])
  surftypes <- c("trend", "cif")[c(trend, cif)]

  # marked point process?
  mrkvals <- attr(x,"mrkvals")
  marked <- (length(mrkvals) > 1)
  if(marked & superimpose) {
    data.types <- levels(data$marks)
    if(any(sort(data.types) != sort(mrkvals)))
      stop(paste("Data marks are different from mark",
                 "values for argument x.\n"))
  }

  # plotting style
  howmat <- outer(how, c("persp", "image", "contour"), "==")
  howmatch <- apply(howmat, 1, any)
  if (any(!howmatch)) 
    stop(paste("unrecognised option", how[!howmatch]))

  # start plotting
  if(pause)
    oldpar <- par(ask = TRUE)
  on.exit(if(pause) par(oldpar))

  
  for(ttt in surftypes) {
    xs <- x[[ttt]]
    for (i in seq(mrkvals)) {
      level <- mrkvals[i]
      main <- paste("Fitted", ttt, "\n",
                    if(marked) paste("mark =", level) else NULL)
      for (style in how) {
        switch(style,
               persp = {
                 do.call("persp",
                         resolve.defaults(list(xs[[i]]),
                                          list(...), 
                                          spatstat.options("par.persp")[[1]],
                                          list(xlab="x", zlab=ttt, main=main)))
               },
               image = {
                 do.call("image",
                         resolve.defaults(list(xs[[i]]),
                                          list(...),
                                          list(main=main)))
                 if(superimpose) {
                   if(marked) plot(data[data$marks == level],
                                   add = TRUE)
                   else plot(data,add=TRUE)
                 }
               },
               contour = {
                 do.call("contour",
                         resolve.defaults(list(xs[[i]]),
                                          list(...),
                                          list(main=main)))
                 if (superimpose) {
                   if(marked) plot(data[data$marks == level],
                                   add = TRUE)
                   else plot(data,add=TRUE)
                 }
               },
               {
                 stop(paste("Unrecognised plot style", style))
               })
      }
    }
  }
  return(invisible())
}

print.plotppm <- function(x, ...) {
  verifyclass(x, "plotppm")
  trend   <- x$trend
  cif     <- x$cif
  mrkvals <- attr(x, "mrkvals")
  ntypes  <- length(mrkvals)
  unmarked <- (ntypes == 1 )
  cat("Object of class \`plotppm\'\n")
  if(unmarked)
    cat("Computed for an unmarked point process\n")
  else {
    cat("Computed for a marked point process, with mark values:\n")
    print(mrkvals)
  }
  cat("Contains the following components:\n")
  if(!is.null(trend)) {
    cat("\n$trend:\tFitted trend.\n")
    if(unmarked) {
      cat("A list containing 1 image\n")
      print(trend[[1]], ...)
    } else {
      cat(paste("A list of", ntypes, "images\n"))
      cat("Typical details:\n")
      print(trend[[1]], ...)
    }
  }
  if(!is.null(cif)) {
    cat("\n$cif:\tFitted conditional intensity.\n")
    if(unmarked) {
      cat("A list containing 1 image\n")
      print(cif[[1]], ...)
    } else {
      cat(paste("A list of", ntypes, "images\n"))
      cat("Typical details:\n")
      print(cif[[1]], ...)
    }
  }
  invisible(NULL)
}

