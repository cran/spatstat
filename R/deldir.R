#
# deldir.R
#
# Interface to deldir package
#
#  $Revision: 1.8 $ $Date: 2011/05/18 01:39:11 $
#

.Spatstat.use.trigraf <- TRUE
.Spatstat.use.trigrafS <- TRUE

dirichlet <- function(X) {
  w <- X$window
  dd <- deldir(X$x, X$y, rw=c(w$xrange,w$yrange))
  pp <- lapply(tile.list(dd), function(z) { owin(poly=z[c("x","y")]) })
  dir <- tess(tiles=pp, window=as.rectangle(w))
  if(w$type != "rectangle")
    dir <- intersect.tess(dir, w)
  return(dir)
}

delaunay <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X)
  nX <- npoints(X)
  if(nX < 3) return(NULL)
  w <- X$window
  dd <- deldir(X$x, X$y, rw=c(w$xrange, w$yrange))
  a <- dd$delsgs[,5]
  b <- dd$delsgs[,6]
  if(.Spatstat.use.trigrafS) {
    # first ensure a[] < b[]
    swap <- (a > b)
    if(any(swap)) {
      oldb <- b
      b[swap] <- a[swap]
      a[swap] <- oldb[swap]
    }
    # next ensure a is sorted
    o <- order(a, b)
    a <- a[o]
    b <- b[o]
    # 
    nv <- nX
    ne <- length(a)
    z <- .C("trigrafS",
            nv = as.integer(nv),
            ne = as.integer(ne),
            ie = as.integer(a - 1),
            je = as.integer(b - 1),
            nt = as.integer(integer(1)),
            it = as.integer(integer(ne)),
            jt = as.integer(integer(ne)),
            kt = as.integer(integer(ne)),
            PACKAGE="spatstat")
    tlist <- with(z, cbind(it, jt, kt)[1:nt, ]) + 1
  } else if(.Spatstat.use.trigraf) {
    nv <- nX
    ne <- length(a)
    z <- .C("trigraf",
            nv = as.integer(nv),
            ne = as.integer(ne),
            ie = as.integer(a - 1),
            je = as.integer(b - 1),
            nt = as.integer(integer(1)),
            it = as.integer(integer(ne)),
            jt = as.integer(integer(ne)),
            kt = as.integer(integer(ne)),
            scratch = as.integer(integer(ne)),
            PACKAGE="spatstat")
    tlist <- with(z, cbind(it, jt, kt)[1:nt, ]) + 1
  } else {
    tlist <- matrix(integer(0), 0, 3)
    for(i in seq_len(nX)) {
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
  }
  # assemble coordinates of triangles
  x <- X$x
  y <- X$y
  xtri <- matrix(x[tlist], nrow(tlist), 3)
  ytri <- matrix(y[tlist], nrow(tlist), 3)
  # ensure triangle vertices are in anticlockwise order
  ztri <- ytri - min(y)
  dx <- cbind(xtri[,2]-xtri[,1], xtri[,3]-xtri[,2], xtri[,1]-xtri[,3])
  zm <- cbind(ztri[,1]+ztri[,2], ztri[,2]+ztri[,3], ztri[,3]+ztri[,1])
  negareas <- apply(dx * zm, 1, sum)
  clockwise <- (negareas > 0)
  #
  if(any(clockwise)) {
    xc <- xtri[clockwise,]
    yc <- ytri[clockwise,]
    xtri[clockwise,] <- xc[,c(1,3,2)]
    ytri[clockwise,] <- yc[,c(1,3,2)]
  }
  # make tile list
  tiles <- list()
  for(m in seq_len(nrow(tlist))) {
    p <- list(x=xtri[m,], y=ytri[m,])
    tiles[[m]] <- owin(poly=p, check=FALSE)
  }

  wc <- convexhull.xy(x, y)
  del <- tess(tiles=tiles, window=wc)
  if(w$type != "rectangle")
    del <- intersect.tess(del, w)
  return(del)
}

  
  
