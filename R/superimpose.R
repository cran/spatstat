# superimpose.R
#
# $Revision: 1.8 $ $Date: 2007/03/16 03:31:37 $
#
# This has been taken out of ppp.S
#
############################# 

"superimpose" <-
  function(..., W=NULL)
{
  # superimpose any number of point patterns
  
  arglist <- list(...)

  if(length(arglist) == 0)
    stop("No point patterns given")
  
  if(length(arglist) == 1 && inherits(arglist[[1]], "list"))
    arglist <- arglist[[1]]

  # determine window
  if(!is.null(W))
    W <- as.owin(W)
  else {
    # extract windows from ppp objects
    isppp <- unlist(lapply(arglist, is.ppp))
    Wlist <- lapply(arglist[isppp], function(x) { x$window })
    # compute bounding boxes of other arguments
    Blist <- lapply(arglist[!isppp], bounding.box.xy)
    Wlist <- append(Wlist, Blist)
    # take the union of all the windows
    W <- Wlist[[1]]
    nW <- length(Wlist)
    if(nW > 1)
      for(i in 2:nW)
        W <- union.owin(W, Wlist[[i]])
  }
     
  # concatenate lists of (x,y) coordinates
  XY <- do.call("concatxy", arglist)

  # create the point pattern
  OUT <- ppp(XY$x, XY$y, window=W)
  
  # find out whether the arguments are marked patterns
  getmarks <- function(x) {
    if(is.ppp(x)) return(marks(x, dfok=FALSE))
    m <- x$marks
    if(is.data.frame(m))
      stop("Sorry, not implemented for data frames of marks")
    return(m)
  }
  Mlist <- lapply(arglist, getmarks)
  ismarked <- !unlist(lapply(Mlist, is.null))
  isfactor <- unlist(lapply(Mlist, is.factor))

  if(any(ismarked) && !all(ismarked))
    warning("Some, but not all, patterns contain marks -- ignored.")
  if(any(isfactor) && !all(isfactor))
    stop("Patterns have incompatible marks - some are factors, some are not")

  if(!all(ismarked)) {
    # Assume all patterns unmarked.
    # If patterns are not named, return the superimposed point pattern.
    nama <- names(arglist)
    if(is.null(nama) || any(nama == ""))
      return(OUT)
    # Patterns are named. Make marks from names.
    len <- unlist(lapply(arglist, function(a) { length(a$x) }))
    M <- factor(rep(nama, len), levels=nama)
    OUT <- OUT %mark% M
    return(OUT)
  }

  # All patterns are marked.
  # Concatenate vectors of marks
  if(!all(isfactor))
    # continuous marks
    M <- unlist(Mlist)
  else {
    # multitype
    Llist <- lapply(Mlist, levels)
    lev <- unique(unlist(Llist))
    codesof <- function(x, lev) { as.integer(factor(x, levels=lev)) }
    Mlist <- lapply(Mlist, codesof, lev=lev)
    M <- factor(unlist(Mlist), levels=codesof(lev,lev), labels=lev)
  }
  OUT <- OUT %mark% M
  return(OUT)
}

