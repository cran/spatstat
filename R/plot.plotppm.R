#
# plot.plotppm.R
#
# engine of plot method for ppm
#
# $Revision: 1.3 $  $Date: 2004/08/30 04:55:22 $
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
      main <- if(marked) paste("mark =", level) else ""
      for (style in how) {
        switch(style,
               persp = {
                 do.call("persp",
                         resolve.defaults(list(xs[[i]]),
                                          list(...), 
                                          spatstat.options("par.persp")[[1]],
                                          list(xlab="x", main=main)))
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
