#
# deldir.R
#
# Interface to deldir package
#
#  $Revision: 1.1 $ $Date: 2008/07/25 23:10:06 $
#

dirichlet <- function(X) {
  w <- X$window
  dd <- deldir(X$x, X$y, rw=c(w$xrange,w$yrange))
  pp <- lapply(tile.list(dd), function(z) { owin(poly=z[c("x","y")]) })
  if(w$type != "rectangle")
    pp <- lapply(pp, intersect.owin, B=w)
  tess(tiles=pp)
}

delaunay <- function(X) {
  stopifnot(is.ppp(X))
  if(X$n < 3) return(NULL)
  w <- X$window
  dd <- deldir(X$x, X$y, rw=c(w$xrange, w$yrange))
  a <- dd$delsgs[,5]
  b <- dd$delsgs[,6]
  tlist <- data.frame(i=integer(0), j=integer(0),k=integer(0))
  for(i in seq(X$n)) {
    # find all Delaunay neighbours of i 
    jj <- c(b[a==i], a[b==i])
    jj <- sort(unique(jj))
    # select those with a higher index than i
    jj <- jj[jj > i]
    # find pairs of neighbours which are Delaunay neighbours
    # (thus, triangles where the first numbered vertex is i)
    if(length(jj) > 0) 
      for(j in jj) {
        kk <- c(b[a == j], a[b == j])
        kk <- kk[(kk %in% jj) & (kk > j)]
        if(length(kk) > 0)
          for(k in kk) 
            # add (i,j,k) to list of triangles (i < j < k)
            tlist <- rbind(tlist, c(i, j, k))
      }
  }
  # make tile list
  tiles <- list()
  for(m in seq(nrow(tlist))) {
    Xijk <- X[as.integer(tlist[m,])]
    p <- list(x=Xijk$x, y=Xijk$y)
    if(area.xypolygon(p) < 0)
      p <- lapply(p, rev)
    tiles[[m]] <- owin(poly=p)
  }
  tess(tiles=tiles)
}

  
  
