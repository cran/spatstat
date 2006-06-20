#
plot.splitppp <- function(x, ..., main, arrange=TRUE) {
  xname <- deparse(substitute(x))
  main <- if(!missing(main)) main else xname
  n <- length(x)
  if(!arrange) {
    lapply(names(x),
           function(l, x, ...){plot(x[[l]], main=l, ...)},
           x=x, ...) 
    return(invisible(NULL))
  }
  nrows <- as.integer(floor(sqrt(n)))
  ncols <- as.integer(ceiling(n/nrows))
  # like mfrow(nrows, ncols) plus a banner at the top
  nblank <- n - ncols * nrows
  mat <- matrix(c(rep(1, ncols), 1 + seq(n), rep(0, nblank)),
                byrow=TRUE, ncol=ncols, nrow=nrows+1)
  layout(mat, heights=c(0.2, rep(1, nrows)))
  opa <- par(mar=rep(0,4))
  plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  text(0,0,main,cex=1.8)
  par(mar=c(2,1,1,2))
  lapply(names(x),
         function(l, x, ...){plot(x[[l]], main=l, ...)},
         x=x, ...)
  # revert
  layout(1)
  par(opa)
  return(invisible(NULL))
}
  
