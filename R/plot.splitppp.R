#
plot.splitppp <- function(x, ..., arrange=TRUE) {
  n <- length(x)
  m <- as.integer(floor(sqrt(n)))
  k <- as.integer(ceiling(n/m))
  if(arrange)
    opa <- par(mfrow=c(k, m))
  lapply(names(x),
         function(l, x, ...){plot(x[[l]], main=l, ...)},
         x=x, ...) 
  if(arrange)
    par(opa)
  return(invisible(NULL))
}
  
