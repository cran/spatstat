#
# lpp.R
#
#  $Revision: 1.11 $   $Date: 2012/04/30 08:49:41 $
#
# Class "lpp" of point patterns on linear networks

lpp <- function(X, L) {
  stopifnot(inherits(L, "linnet"))
  if(!is.ppp(X))
    X <- as.ppp(X, W=L$window)
  X <- project2segment(X, as.psp(L))$Xproj
  out <- ppx(data=as.data.frame(X),
             domain=L,
             spatial=1:2)
  class(out) <- c("lpp", class(out))
  return(out)
}

print.lpp <- function(x, ...) {
  stopifnot(inherits(x, "lpp"))
  cat("Point pattern on linear network\n")
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
  print(x$domain, ...)
  invisible(NULL)
}

# plot.lpp removed: plot.ppx sufficient

summary.lpp <- function(object, ...) {
  stopifnot(inherits(object, "lpp"))
  L <- object$domain
  npoints <- nrow(object$data)
  totlen <-  sum(lengths.psp(L$lines))
  out <- list(npoints=npoints,
              totlength=totlen,
              intensity=npoints/totlen,
              nvert=L$vertices$n,
              nedge=L$lines$n,
              unitinfo=summary(unitname(L)))
  class(out) <- "summary.lpp"
  return(out)
}

print.summary.lpp <- function(x, ...) {
  cat("Point pattern on linear network\n")
  cat(paste(x$npoints, "points\n"))
  cat(paste("Linear network with",
            x$nvert, "vertices and",
            x$nedge, "edges\n"))
  u <- x$unitinfo
  cat(paste("Total edge length", x$totlength, u$plural, u$explain, "\n"))
  cat(paste("Average intensity", x$intensity,
            "points per", if(u$vanilla) "unit length" else u$singular, "\n"))
  invisible(NULL)
}

as.ppp.lpp <- function(X, ..., fatal=TRUE) {
  verifyclass(X, "lpp", fatal=fatal)
  L <- X$domain
  as.ppp(as.data.frame(X), W=L$window)
}

as.linnet.lpp <- function(X, ..., fatal=TRUE) {
  verifyclass(X, "lpp", fatal=fatal)
  X$domain
}
  
"[.lpp" <- function (x, i, ...) {
  # invoke [.ppx
  y <- NextMethod("[")
  class(y) <- c("lpp", class(y))
  return(y)
}

unitname.lpp <- function(x) {
  u <- unitname(x$domain)
  return(u)
}

"unitname<-.lpp" <- function(x, value) {
  w <- x$domain
  unitname(w) <- value
  x$domain <- w
  return(x)
}

"marks<-.lpp" <- function(x, ..., value) {
  Y <- NextMethod("marks<-")
  class(Y) <- c("lpp", class(Y))
  Y
}
  
unmark.lpp <- function(X) {
  Y <- NextMethod("unmark")
  class(Y) <- c("lpp", class(Y))
  Y
}

as.psp.lpp <- function(x, ..., fatal=TRUE){
  verifyclass(x, "lpp", fatal=fatal)
  return(x$domain$lines)
}

