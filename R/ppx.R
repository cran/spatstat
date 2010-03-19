#
#   ppx.R
#
#  class of general point patterns in any dimension
#
#  $Revision: 1.4 $  $Date: 2010/03/18 02:45:45 $
#

ppx <- function(data, domain=NULL, spatial=NULL, temporal=NULL) {
  data <- as.hyperframe(data)
  suitable <- (summary(data)$storage == "dfcolumn")
  if(is.null(spatial) && is.null(temporal)) {
    # assume all suitable columns of data are spatial coordinates
    spatial <- suitable
    temporal <- FALSE & suitable
  } else {
    # 'spatial' determines which columns contain spatial coordinates
    # 'temporal' determines which columns contain temporal coordinates
    # convert them to logical vectors
    vnames <- names(data)
    names(vnames) <- vnames
    spatial <- vnames %in% vnames[spatial]
    temporal <- vnames %in% vnames[temporal]
    # check for conflicts
    if(any(nbg <- spatial & temporal)) 
      stop(paste(ngettext(sum(nbg), "coordinate", "coordinates"),
                 commasep(sQuote(vnames[spatial & temporal])),
                 "selected as both spatial and temporal"))
    # check that the coordinate columns contain ordinary data
    if(any(nbg <- (spatial & !suitable)))
      stop(paste(ngettext(sum(nbg), "column", "columns"),
                 commasep(dQuote(vnames[nbg])),
                 ngettext(sum(nbg), "contain", "contains"),
                 "data that are not suitable as spatial coordinates"))
    if(any(nbg <- (temporal & !suitable)))
      stop(paste(ngettext(sum(nbg), "column", "columns"),
                 commasep(dQuote(vnames[nbg])),
                 ngettext(sum(nbg), "contain", "contains"),
                 "data that are not suitable as time coordinates"))
  }
  ctype <- ifelse(spatial, "spatial", ifelse(temporal, "temporal", "mark"))
  ctype <- factor(ctype, levels=c("spatial", "temporal", "mark"))
  out <- list(data=data, ctype=ctype, domain=domain)
  class(out) <- "ppx"
  return(out)
}

is.ppx <- function(x) { inherits(x, "ppx") }

npoints.ppx <- function(x) { nrow(x$data) }

print.ppx <- function(x, ...) {
  cat("Multidimensional point pattern\n")
  sd <- summary(x$data)
  np <- sd$ncases
  nama <- sd$col.names
  cat(paste(np, ngettext(np, "point", "points"), "\n"))
  if(any(iscoord <- (x$ctype == "spatial")))
    cat(paste(sum(iscoord), "-dimensional space coordinates ",
              paren(paste(nama[iscoord], collapse=",")), "\n", sep=""))
  if(any(istime <- (x$ctype == "temporal")))
    cat(paste(sum(istime), "-dimensional time coordinates ",
              paren(paste(nama[istime], collapse=",")), "\n", sep=""))
  if(any(ismark <- (x$ctype == "mark"))) 
    cat(paste(sum(ismark), ngettext(sum(ismark), "column", "columns"),
              "of marks:",
              commasep(sQuote(nama[ismark])), "\n"))
  if(!is.null(x$domain)) {
    cat("Domain:\n\t")
    print(x$domain)
  }
  invisible(NULL)
}

plot.ppx <- function(x, ...) {
  xname <- deparse(substitute(x))
  coo <- coords(x)
  if(ncol(coo) != 2)
    stop(paste("Don't know how to plot a general point pattern in",
               ncol(coo), "dimensions"))
  do.call("plot", resolve.defaults(list(x$domain),
                                   list(...),
                                   list(main=xname)))
  do.call.matched("points",
                  append(list(x=coo), list(...)),
                  extrargs=c("type", "pch", "col", "bg", "cex", "lty", "lwd"))
  return(invisible(NULL))
}

coords <- function(x, ...) {
  UseMethod("coords")
}

coords.ppx <- function(x, ..., spatial=TRUE, temporal=TRUE) {
  ctype <- x$ctype
  chosen <- (ctype == "spatial" & spatial) | (ctype == "temporal" & temporal)
  as.data.frame(x$data[, chosen])
}

coords.ppp <- function(x, ...) { data.frame(x=x$x,y=x$y) }

as.hyperframe.ppx <- function(x, ...) { x$data }

as.data.frame.ppx <- function(x, ...) { as.data.frame(x$data) } 

marks.ppx <- function(x, ..., drop=TRUE) {
  ctype <- x$ctype
  chosen <- (ctype == "mark")
  x$data[, chosen, drop=drop]
}

"marks<-.ppx" <- function(x, ..., value) {
  ctype <- x$ctype
  retain <- (ctype != "mark")
  coorddata <- x$data[, retain, drop=TRUE]
  if(!is.data.frame(value) && !is.hyperframe(value))
    value <- data.frame(marks=value)
  newdata <- cbind(coorddata, value)
  newctype <- c(ctype, rep("mark", ncol(value)))
  out <- list(data=newdata, ctype=newctype, domain=x$domain)
  class(out) <- "ppx"
  return(out)
}

  

