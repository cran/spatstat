#
# summary.quad.R
#
#  summary() method for class "quad"
#
#  $Revision: 1.3 $ $Date: 2004/01/27 07:06:14 $
#
summary.quad <- function(object, ...) {
  verifyclass(object, "quad")
  s <- list(
       data  = summary.ppp(object$data),
       dummy = summary.ppp(object$dummy),
       param = object$param)
  doit <- function(ww) {
    return(list(range=range(ww), sum=sum(ww)))
  }
  w <- object$w
  Z <- is.data(object)
  s$w <- list(all=doit(w), data=doit(w[Z]), dummy=doit(w[!Z]))
  class(s) <- "summary.quad"
  return(s)
}

print.summary.quad <- function(x, ..., dp=3) {
  cat("Quadrature scheme = data + dummy + weights\n")
  pa <- x$param
  if(is.null(pa))
    cat("created by an unknown function.\n")
  cat("Data pattern:\n")
  print(x$data, dp=dp)

  cat("\n\nDummy quadrature points:\n")
  # How they were computed
  if(!is.null(pa)) {
    dumpar <- pa$dummy
    if(is.null(dumpar))
      cat("(provided manually)\n")
    else if(!is.null(dumpar$nd)) 
      cat(paste("(", dumpar$nd[1], "x", dumpar$nd[2],
                "grid, plus 4 corner points)\n"))
    else
      cat("(rule for creating dummy points not understood)")
  }
  # Description of them
  print(x$dummy, dp=dp)

  cat("\n\nQuadrature weights:\n")
  # How they were computed
  if(!is.null(pa)) {
    wpar <- pa$weight
    if(is.null(wpar))
      cat("(values provided manually)\n")
    else if(!is.null(wpar$method)) {
      if(wpar$method=="grid") {
        cat(paste("(counting weights based on",
                  wpar$ntile[1], "x", wpar$ntile[2],
                  "array of rectangular tiles)\n"))
      } else if(wpar$method=="dirichlet") {
        cat(paste("(Dirichlet tile areas, computed",
                  if(wpar$exact) "exactly" else "by pixel approximation",
                  ")\n"))
      } else
      cat("(rule for creating dummy points not understood)\n")
    }
  }
  # Description of them
  doit <- function(ww) {
    cat(paste("range: ",
              "[",
              paste(signif(ww$range, digits=dp), collapse=", "),
              "]\t",
              "total: ",
              signif(ww$sum, digits=dp),
              "\n", sep=""))
  }
  cat("All weights:\n\t")
  doit(x$w$all)
  cat("Weights on data points:\n\t")
  doit(x$w$data)
  cat("Weights on dummy points:\n\t")
  doit(x$w$dummy)

  return(invisible(NULL))
}

    
print.quad <- function(x, ...) {
  cat("Quadrature scheme\n")
  cat(paste(x$data$n, "data points, ", x$dummy$n, "dummy points\n"))
  cat(paste("Total weight ", sum(x$w), "\n"))
  return(invisible(NULL))
}
