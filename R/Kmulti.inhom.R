#
#	Kmulti.inhom.S		
#
#	$Revision: 1.22 $	$Date: 2009/06/10 01:26:28 $
#
#
# ------------------------------------------------------------------------

"Kcross.inhom" <- 
function(X, i, j, lambdaI, lambdaJ, ...,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
         lambdaIJ=NULL)
{
  verifyclass(X, "ppp")
  if(!is.marked(X))
    stop("point pattern has no marks (no component 'marks')")
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  J <- (marx == j)
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points with mark j =", j)
  result <- Kmulti.inhom(X, I, J, lambdaI, lambdaJ, ...,
                         r=r,breaks=breaks,correction=correction,
                         lambdaIJ=lambdaIJ, Iname=Iname, Jname=Jname)
  result <-
    rebadge.fv(result,
               substitute(Kinhom[i,j](r), list(i=paste(i), j=paste(j))),
               "K*i",
               new.yexp=substitute(Kinhom[list(i,j)](r),
                                   list(i=paste(i), j=paste(j))))
  return(result)
}

"Kdot.inhom" <- 
function(X, i, lambdaI, lambdadot, ...,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
         lambdaIdot=NULL)
{
  verifyclass(X, "ppp")
  if(!is.marked(X))
    stop("point pattern has no marks (no component 'marks')")
  if(missing(correction))
    correction <- NULL

  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]

  I <- (marx == i)
  J <- rep(TRUE, X$n)  # i.e. all points
  Iname <- paste("points with mark i =", i)
  Jname <- paste("points")
	
  result <- Kmulti.inhom(X, I, J, lambdaI, lambdadot, ...,
                         r=r,breaks=breaks,correction=correction,
                         lambdaIJ=lambdaIdot,
                         Iname=Iname, Jname=Jname)
  result <-
    rebadge.fv(result,
               substitute(Kdot.inhom[i](r), list(i=paste(i))),
               "K.i")
  return(result)
}


"Kmulti.inhom"<-
function(X, I, J, lambdaI, lambdaJ, 
         ...,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
         lambdaIJ=NULL,
         Iname = "points satisfying condition I",
         Jname = "points satisfying condition J")
{
  verifyclass(X, "ppp")

  extras <- list(...)
  if(length(extras) > 0)
    warning(paste("Unrecognised arguments", names(extras)))
        
  npoints <- X$n
  x <- X$x
  y <- X$y
  W <- X$window
  area <- area.owin(W)

  # validate edge correction
  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("border", "isotropic", "Ripley", "translate")

  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  # validate I, J 
  if(!is.logical(I) || !is.logical(J))
    stop("I and J must be logical vectors")
  if(length(I) != npoints || length(J) != npoints)
    stop(paste("The length of I and J must equal",
               "the number of points in the pattern"))
	
  nI <- sum(I)
  nJ <- sum(J)
  if(nI == 0) stop(paste("There are no", Iname))
  if(nJ == 0) stop(paste("There are no", Jname))

  # r values 
  rmaxdefault <- rmax.rule("K", W, nJ/area)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  # intensity data
  if(is.im(lambdaI)) {
    # look up intensity values
    lambdaI <- safelookup(lambdaI, X[I])
  } else if(is.vector(lambdaI) && is.numeric(lambdaI)) {
    # validate intensity vector
    if(length(lambdaI) != nI)
      stop(paste("The length of", sQuote("lambdaI"),
                 "should equal the number of", Iname))
  } else 
  stop(paste(sQuote("lambdaI"), "should be a vector or an image"))

  if(is.im(lambdaJ)) {
    # look up intensity values
    lambdaJ <- safelookup(lambdaJ, X[J])
  } else if(is.vector(lambdaJ) && is.numeric(lambdaJ)) {
    # validate intensity vector
    if(length(lambdaJ) != nJ)
      stop(paste("The length of", sQuote("lambdaJ"),
                 "should equal the number of", Jname))
  } else 
  stop(paste(sQuote("lambdaJ"), "should be a vector or an image"))

  # Weight for each pair
  if(!is.null(lambdaIJ)) {
    if(!is.matrix(lambdaIJ))
      stop("lambdaIJ should be a matrix")
    if(nrow(lambdaIJ) != nI)
      stop(paste("nrow(lambdaIJ) should equal the number of", Iname))
    if(ncol(lambdaIJ) != nJ)
      stop(paste("ncol(lambdaIJ) should equal the number of", Jname))
  }

  # Recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))
        
  # this will be the output data frame
  # It will be given more columns later
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", substitute(Kmulti.inhom(r), NULL),
          "theo", , alim, c("r","%spois(r)"), desc, fname="Kmi")

# identify close pairs of points
  XI <- X[I]
  XJ <- X[J]
  close <- crosspairs(XI, XJ, max(r))
# map (i,j) to original serial numbers in X
  orig <- seq(npoints)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
# eliminate any identical pairs
  if(any(I & J)) {
    ok <- (iX != jX)
    if(!all(ok)) {
      close$i  <- close$i[ok]
      close$j  <- close$j[ok]
      close$xi <- close$xi[ok]
      close$yi <- close$yi[ok]
      close$xj <- close$xj[ok]
      close$yj <- close$yj[ok]
      close$dx <- close$dx[ok]
      close$dy <- close$dy[ok]
      close$d  <- close$d[ok]
    }
  }
# extract information for these pairs (relative to orderings of XI, XJ)
  dclose <- close$d
  icloseI  <- close$i
  jcloseJ  <- close$j
        
# Form weight for each pair
  if(is.null(lambdaIJ))
    weight <- 1/(lambdaI[icloseI] * lambdaJ[jcloseJ])
  else 
    weight <- 1/lambdaIJ[cbind(icloseI, jcloseJ)]

# Compute estimates by each of the selected edge corrections.

  if(any(correction == "border" | correction == "bord.modif")) {
    # border method
    # Compute distances to boundary
    b <- bdist.points(XI)
    bI <- b[icloseI]
    # apply reduced sample algorithm
    RS <- Kwtsum(dclose, bI, weight, b, 1/lambdaI, breaks)
    if(any(correction == "border")) {
      Kb <- RS$ratio
      K <- bind.fv(K, data.frame(border=Kb), "%sbord(r)",
                   "border-corrected estimate of %s",
                   "border")
    }
    if(any(correction == "bord.modif")) {
      Kbm <- RS$numerator/eroded.areas(W, r)
            K <- bind.fv(K, data.frame(bord.modif=Kbm), "%sbord*(r)",
                         "modified border-corrected estimate of %s",
                         "bord.modif")
          }
        }
        if(any(correction == "translate")) {
          # translation correction
            edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
            allweight <- edgewt * weight
            wh <- whist(dclose, breaks$val, allweight)
            Ktrans <- cumsum(wh)/area
            rmax <- diameter(W)/2
            Ktrans[r >= rmax] <- NA
            K <- bind.fv(K, data.frame(trans=Ktrans), "%strans(r)",
                         "translation-corrected estimate of %s",
                         "trans")
        }
        if(any(correction == "isotropic")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
            allweight <- edgewt * weight
            wh <- whist(dclose, breaks$val, allweight)
            Kiso <- cumsum(wh)/area
            rmax <- diameter(W)/2
            Kiso[r >= rmax] <- NA
            K <- bind.fv(K, data.frame(iso=Kiso), "%siso(r)",
                         "Ripley isotropic correction estimate of %s",
                         "iso")
        }
        # default is to display them all
        attr(K, "fmla") <- . ~ r
        unitname(K) <- unitname(X)
        return(K)
          
}
