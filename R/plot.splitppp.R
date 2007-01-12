#
plot.listof <- plot.splitppp <- function(x, ..., main, arrange=TRUE,
                                         nrows=NULL, ncols=NULL) {
  xname <- deparse(substitute(x))

  if(is.null(names(x)))
    names(x) <- paste("Component_", seq(length(x)), sep="")

  if(!arrange) {
    lapply(names(x),
           function(l, x, ...){plot(x[[l]], main=l, ...)},
           x=x, ...) 
    return(invisible(NULL))
  }

  # decide whether to plot a main header
  main <- if(!missing(main)) main else xname
  banner <- (sum(nchar(as.character(main))) > 0)
  if(length(main) > 1)
    main <- paste(main, collapse="\n")
  nlines <- if(!is.character(main)) 1 else length(unlist(strsplit(main, "\n")))
  # determine arrangement of plots
  # arrange like mfrow(nrows, ncols) plus a banner at the top
  n <- length(x)
  if(is.null(nrows) && is.null(ncols)) {
    nrows <- as.integer(floor(sqrt(n)))
    ncols <- as.integer(ceiling(n/nrows))
  } else if(!is.null(nrows) && is.null(ncols))
    ncols <- as.integer(ceiling(n/nrows))
  else if(is.null(nrows) && !is.null(ncols))
    nrows <- as.integer(ceiling(n/ncols))
  else stopifnot(nrows * ncols >= length(x))
  nblank <- ncols * nrows - n
  # declare layout
  mat <- matrix(c(seq(n), rep(0, nblank)),
                byrow=TRUE, ncol=ncols, nrow=nrows)
  heights <- rep(1, nrows)
  if(banner) {
    # Increment existing panel numbers
    # New panel 1 is the banner
    panels <- (mat > 0)
    mat[panels] <- mat[panels] + 1
    mat <- rbind(rep(1,ncols), mat)
    heights <- c(0.1 * (1 + nlines), heights)
  }
  layout(mat, heights=heights)
  # plot banner
  if(banner) {
    opa <- par(mar=rep(0,4), xpd=TRUE)
    plot(numeric(0),numeric(0),type="n",ann=FALSE,axes=FALSE,
         xlim=c(-1,1),ylim=c(-1,1))
    cex <- resolve.defaults(list(...), list(cex.title=2))$cex.title
    text(0,0,main, cex=cex)
  }
  # plot panels
  npa <- par(mar=c(2,1,1,2))
  if(!banner) opa <- npa
  lapply(names(x),
         function(l, x, ...){plot(x[[l]], main=l, ...)},
         x=x, ...)
  # revert
  layout(1)
  par(opa)
  return(invisible(NULL))
}
  
density.splitppp <- function(x, ...) {
  u <- lapply(x, density, ...)
  class(u) <- c("listof", class(u))
  return(u)
}
