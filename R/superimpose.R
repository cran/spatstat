# superimpose.R
#
# $Revision: 1.12 $ $Date: 2008/10/15 00:16:35 $
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
    isnull <- unlist(lapply(arglist, is.null))
    if(any(!isppp & !isnull)) {
      Blist <- lapply(arglist[!isppp & !isnull], bounding.box.xy)
      Bisnull <- unlist(lapply(Blist, is.null))
      Wlist <- append(Wlist, Blist[!Bisnull])
    }
    # take the union of all the windows
    W <- Wlist[[1]]
    nW <- length(Wlist)
    if(nW > 1)
      for(i in 2:nW)
        W <- union.owin(W, Wlist[[i]])
  }
     
  # concatenate lists of (x,y) coordinates
  XY <- do.call("concatxy", arglist)

  # create the point pattern without checking
  OUT <- ppp(XY$x, XY$y, window=W, check=FALSE)
  
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
      return(as.ppp(OUT, check=TRUE))
    # Patterns are named. Make marks from names.
    len <- unlist(lapply(arglist, function(a) { length(a$x) }))
    M <- factor(rep(nama, len), levels=nama)
    OUT <- OUT %mark% M
    return(as.ppp(OUT, check=TRUE))
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
  return(as.ppp(OUT, check=TRUE))
}

