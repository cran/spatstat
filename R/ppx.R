#
#   ppx.R
#
#  class of general point patterns in any dimension
#
#  $Revision: 1.21 $  $Date: 2010/08/17 07:20:43 $
#

ppx <- function(data, domain=NULL, spatial=NULL, temporal=NULL) {
  data <- as.hyperframe(data)
  suitable <- (summary(data)$storage == "dfcolumn")
  if(is.null(spatial) && is.null(temporal)) {
    # assume all suitable columns of data are spatial coordinates
    isspatial <- suitable
    istemporal <- FALSE & suitable
  } else {
    # 'spatial' determines which columns contain spatial coordinates
    # 'temporal' determines which columns contain temporal coordinates
    # convert them to logical vectors
    vnames <- names(data)
    names(vnames) <- vnames
    isspatial <- vnames %in% vnames[spatial]
    istemporal <- vnames %in% vnames[temporal]
    if(is.null(spatial))
      isspatial <- suitable & !istemporal
    # check for conflicts
    if(any(nbg <- isspatial & istemporal)) 
      stop(paste(ngettext(sum(nbg), "coordinate", "coordinates"),
                 commasep(sQuote(vnames[isspatial & istemporal])),
                 "selected as both spatial and temporal"))
    # check that the coordinate columns contain ordinary data
    if(any(nbg <- (isspatial & !suitable)))
      stop(paste(ngettext(sum(nbg), "column", "columns"),
                 commasep(dQuote(vnames[nbg])),
                 ngettext(sum(nbg), "contain", "contains"),
                 "data that are not suitable as spatial coordinates"))
    if(any(nbg <- (istemporal & !suitable)))
      stop(paste(ngettext(sum(nbg), "column", "columns"),
                 commasep(dQuote(vnames[nbg])),
                 ngettext(sum(nbg), "contain", "contains"),
                 "data that are not suitable as time coordinates"))
  }
  ctype <- ifelse(isspatial, "spatial", ifelse(istemporal, "temporal", "mark"))
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
  coo <- coords(x, temporal=FALSE)
  if(any(istime <- (x$ctype == "temporal")))
    warning(paste(sum(istime), "time coordinate(s) ignored"))
  dom <- x$domain
  m <- ncol(coo)
  if(m == 1) {
    coo <- as.vector(coo)
    ran <- diff(range(coo))
    ylim <- c(-1,1) * ran/20
    do.call("plot.default",
            resolve.defaults(list(coo, rep(0, length(coo))),
                             list(...),
                             list(asp=1, ylim=ylim,
                                  axes=FALSE, xlab="", ylab="")))
    axis(1, pos=ylim[1])
  } else if(m == 2) {
    if(is.null(dom)) {
      # plot x, y coordinates only
      do.call.matched("plot.default",
                      resolve.defaults(list(coo[,1], coo[,2], asp=1),
                                       list(...),
                                       list(main=xname)))
    } else {
      # plot domain, whatever it is
      do.call("plot", resolve.defaults(list(dom),
                                       list(...),
                                       list(main=xname)))
      # convert to ppp
      x2 <- ppp(coo[,1], coo[,2], window=as.owin(dom),
                marks=as.data.frame(marks(x)), check=FALSE)
      # invoke plot.ppp
      return(do.call("plot", resolve.defaults(list(x2),
                                              list(add=TRUE),
                                              list(...))))
    }
  } else if(m == 3) {
    # convert to pp3
    x3 <- pp3(coo[,1], coo[,2], coo[,3], domain=dom, check=FALSE)
    # invoke plot.pp3
    do.call("plot", resolve.defaults(list(x3),
                                     list(...),
                                     list(main=xname)))
  } else stop(paste("Don't know how to plot a general point pattern in",
               ncol(coo), "dimensions"))
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

as.data.frame.ppx <- function(x, ...) { as.data.frame(x$data, ...) } 

as.matrix.ppx <- function(x, ...) { as.matrix(as.data.frame(x, ...)) }

marks.ppx <- function(x, ..., drop=TRUE) {
  ctype <- x$ctype
  chosen <- (ctype == "mark")
  x$data[, chosen, drop=drop]
}

"marks<-.ppx" <- function(x, ..., value) {
  ctype <- x$ctype
  retain <- (ctype != "mark")
  coorddata <- x$data[, retain, drop=TRUE]
  if(is.null(value)) {
    newdata <- coorddata
    newctype <- ctype[retain]
  } else {
    if(!is.data.frame(value) && !is.hyperframe(value))
      value <- data.frame(marks=value)
    if(ncol(value) == 0) {
      newdata <- coorddata
      newctype <- ctype[retain]
    } else {
      newdata <- cbind(coorddata, value)
      newctype <- c(ctype, rep("mark", ncol(value)))
    }
  }
  out <- list(data=newdata, ctype=newctype, domain=x$domain)
  class(out) <- "ppx"
  return(out)
}

unmark.ppx <- function(X) {
  marks(X) <- NULL
  return(X)
}

boxx <- function(..., unitname=NULL) {
  if(length(list(...)) == 0)
    stop("No data")
  ranges <- data.frame(...)
  nama <- names(list(...))
  if(is.null(nama) || !all(nzchar(nama)))
    names(ranges) <- paste("x", 1:ncol(ranges),sep="")
  if(nrow(ranges) != 2)
    stop("Data should be vectors of length 2")
  if(any(unlist(lapply(ranges, diff)) <= 0))
    stop("Illegal range: Second element <= first element")
  out <- list(ranges=ranges, units=as.units(unitname))
  class(out) <- "boxx"
  return(out)
}

print.boxx <- function(x, ...) {
  m <- ncol(x$ranges)
  cat(paste(m, "-dimensional box:\n", sep=""))
  bracket <- function(z) paste("[",
                               paste(signif(z, 5), collapse=", "),
                               "]", sep="")
  v <- paste(unlist(lapply(x$ranges, bracket)), collapse=" x ")
  s <- summary(unitname(x))
  cat(paste(v, s$plural, s$explain, "\n"))
  invisible(NULL)
}

unitname.boxx <- function(x) { x$units }

"unitname<-.boxx" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

unitname.ppx <- function(x) { unitname(x$domain) }

"unitname<-.ppx" <- function(x, value) {
  d <- x$domain
  unitname(d) <- value
  x$domain <- d
  return(x)
}

volume.boxx <- function(x) {
  stopifnot(inherits(x, "boxx"))
  prod(unlist(lapply(x$ranges, diff)))
}

diameter.boxx <- function(x) {
  stopifnot(inherits(x, "boxx"))
  sqrt(sum(unlist(lapply(x$ranges, diff))^2))
}

shortside.boxx <- function(x) {
  stopifnot(inherits(x, "boxx"))
  min(unlist(lapply(x$ranges, diff)))
}

eroded.volumes.boxx <- function(x, r) {
  stopifnot(inherits(x, "boxx"))
  ero <- sapply(x$ranges, function(z, r) { pmax(0, diff(z) - 2 * r)}, r=r)
  apply(ero, 1, prod)
}

runifpointx <- function(n, domain) {
  stopifnot(inherits(domain, "boxx"))
  coo <- lapply(domain$ranges,
                function(ra, n) { runif(n, min=ra[1], max=ra[2]) },
                n=n)
  df <- do.call("data.frame", coo)
  ppx(df, domain)
}

rpoisppx <- function(lambda, domain) {
  stopifnot(inherits(domain, "boxx"))
  vol <- volume.boxx(domain)
  stopifnot(is.numeric(lambda) && length(lambda) == 1 && lambda >= 0)
  n <- rpois(1, lambda * vol)
  runifpointx(n, domain)
}
