#
plot.listof <- plot.splitppp <- function(x, ..., main, arrange=TRUE,
                                         nrows=NULL, ncols=NULL,
                                         main.panel=NULL,
                                         mar.panel=c(2,1,1,2)) {
  xname <- deparse(substitute(x))
  n <- length(x)

  if(is.null(names(x)))
    names(x) <- paste("Component_", 1:n, sep="")

  if(is.null(main.panel))
    main.panel <- names(x)
  else {
    stopifnot(is.character(main.panel))
    nmp <- length(main.panel)
    if(nmp == 1)
      main.panel <- rep(main.panel, n)
    else if(nmp != n)
      stop("Incorrect length for main.panel")
  }
  
  if(!arrange) {
    for(i in 1:n) 
      plot(x[[i]], main=main.panel[i], ...)
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
  npa <- par(mar=mar.panel)
  if(!banner) opa <- npa
  for(i in 1:n) 
    plot(x[[i]], main=main.panel[i], ...)
  # revert
  layout(1)
  par(opa)
  return(invisible(NULL))
}
  
density.splitppp <- function(x, ...) {
  as.listof(lapply(x, density, ...))
}
