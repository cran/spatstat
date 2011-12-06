#
#   edges2triangles.R
#
#   $Revision: 1.7 $  $Date: 2011/11/12 10:28:54 $
#

edges2triangles <- function(iedge, jedge, nvert=max(iedge, jedge),
                            ..., check=TRUE, friendly=rep(TRUE, nvert)) {
  usefriends <- !missing(friendly)
  if(check) {
    stopifnot(length(iedge) == length(jedge))
    stopifnot(all(iedge > 0))
    stopifnot(all(jedge > 0))
    if(!missing(nvert)) {
      stopifnot(all(iedge <= nvert))
      stopifnot(all(jedge <= nvert))
    }
    if(usefriends) {
      stopifnot(is.logical(friendly))
      stopifnot(length(friendly) == nvert)
      usefriends <- !all(friendly)
    }
  }
  # zero length data
  if(length(iedge) == 0) return(matrix(, nrow=0, ncol=3))
  # sort in increasing order of 'iedge'
  oi <- order(iedge)
  iedge <- iedge[oi]
  jedge <- jedge[oi]
  # convert to C indexing
  ii <- iedge - 1
  jj <- jedge - 1
  # call C
  storage.mode(nvert) <- storage.mode(ii) <- storage.mode(jj) <- "integer"
  if(!usefriends) {
    zz <- .Call("triograph", nv=nvert, iedge=ii, jedge=jj, PACKAGE="spatstat")
  } else {
    fr <- as.logical(friendly)
    storage.mode(fr) <- "integer"
    zz <- .Call("trioxgraph", nv=nvert, iedge=ii, jedge=jj, friendly=fr,
                PACKAGE="spatstat")
  }
  # convert back to R indexing
  mat <- as.matrix(as.data.frame(zz)) + 1
  return(mat)
}

# compute triangle diameters as well

trianglediameters <- function(iedge, jedge, edgelength, ..., 
                              nvert=max(iedge, jedge), check=TRUE) {
  if(check) {
    stopifnot(length(iedge) == length(jedge))
    stopifnot(length(iedge) == length(edgelength))
    stopifnot(all(iedge > 0))
    stopifnot(all(jedge > 0))
    if(!missing(nvert)) {
      stopifnot(all(iedge <= nvert))
      stopifnot(all(jedge <= nvert))
    }
  }
  # zero length data
  if(length(iedge) == 0)
    return(data.frame(i=integer(0),
                      j=integer(0),
                      k=integer(0),
                      diam=numeric(0)))
  # convert to C indexing
  ii <- iedge - 1
  jj <- jedge - 1
  eij <- edgelength
  # call C
  storage.mode(nvert) <- storage.mode(ii) <- storage.mode(jj) <- "integer"
  storage.mode(eij) <- "double"
  zz <- .Call("triDgraph", nv=nvert, iedge=ii, jedge=jj, edgelength=eij,
              PACKAGE="spatstat")
  df <- as.data.frame(zz)
  colnames(df) <- c("i", "j", "k", "diam")
  # convert back to R indexing
  df[, 1:3] <- df[, 1:3] + 1
  return(df)
}
