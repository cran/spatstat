#
#
#   rmhstart.R
#
#   $Revision: 1.4 $  $Date: 2005/03/10 17:29:47 $
#
#

rmhstart <- function(start, ...) {
  UseMethod("rmhstart")
}

rmhstart.rmhstart <- function(start, ...) {
  return(start)
}

rmhstart.list <- function(start, ...) {
  argnames <- c("x.start","n.start","iseed")
  ok <- argnames %in% names(start) 
  st <- do.call("rmhstart.default", start[argnames[ok]])
  return(st)
}

rmhstart.default <- function(start=NULL, ..., n.start=NULL, x.start=NULL,
                             iseed=NULL)
{
 if(!is.null(start) || length(list(...)) > 0)
    stop("Syntax should be rmhstart(n.start, x.start, iseed)")
 
  ngiven <- !is.null(n.start)
  xgiven <- !is.null(x.start)
  
  # n.start and x.start are incompatible
  if(ngiven && xgiven)
    stop("Give only one of the arguments n.start and x.start")
  if(!ngiven && !xgiven)
    stop("No starting state specified.")

  given <- if(ngiven) "n" else if(xgiven) "x" else "none"

  # Validate arguments
  if(ngiven && !is.numeric(n.start))
      stop("n.start should be numeric")
  if(xgiven) {
    # We can't check x.start properly because we don't have the relevant window
    # Just check that it is INTERPRETABLE as a point pattern  
    xx <- as.ppp(x.start, W=tightbox.xy, fatal=FALSE)
    if(is.null(xx))
      stop(paste("x.start should be a point pattern object,",
                 "or coordinate data in a format recognised by as.ppp"))
  } else
     xx <- NULL

 seed <- if(is.null(iseed)) NULL else rmhseed(iseed)
 
###################################################################
# return augmented list  
  out <- list(n.start=n.start,
              x.start=x.start,
              seed=seed,
              given=given,
              xx=xx)
  class(out) <- c("rmhstart", class(out))
  return(out)
  }

print.rmhstart <- function(x, ...) {
  verifyclass(x, "rmhstart")

  cat("Metropolis-Hastings algorithm starting parameters\n")
  cat("Initial state: ")
  switch(x$given,
         none={ cat("not given\n") },
         x = {
               cat("given as x.start\n")
               if(is.ppp(x$x.start)) 
                 print(x$x.start)
               else
                 cat(paste("(x,y) coordinates of", x$xx$n,
                           "points (window unspecified)\n"))
               cat("\n")
             },
         n = {
           n.start <- x$n.start
           nstring <-
             if(length(n.start) == 1)
               paste(n.start)
             else 
               paste("(", paste(n.start, collapse=","), ")", sep="")
           cat(paste("number fixed at n.start =", nstring, "\n")) }
         )
  if(is.null(x$seed))
    cat("Random number seeds: not fixed (random state)\n")
  else
    print(x$seed)
}

