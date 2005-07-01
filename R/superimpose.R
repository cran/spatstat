# superimpose.R
#
# $Revision: 1.5 $ $Date: 2005/07/01 07:30:25 $
#
# This has been taken out of ppp.S
#
############################# 

"superimpose" <-
  function(...)
{
  # superimpose any number of point patterns
  # ASSUMED TO BE IN THE SAME WINDOW
  
  arglist <- list(...)

  if(length(arglist) == 1 && inherits(arglist[[1]], "list"))
    arglist <- arglist[[1]]
  
  # concatenate lists of (x,y) coordinates
  XY <- do.call("concatxy", arglist)

  # determine window
  P <- arglist[[1]]
  if(!verifyclass(P, "ppp", fatal=FALSE))
    stop("The first argument is not a point pattern object")
  OUT <- ppp(XY$x, XY$y, window=P$window)
  
  # find out whether the arguments are marked patterns
  Mlist <- lapply(arglist, function(x) {x$marks})
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

