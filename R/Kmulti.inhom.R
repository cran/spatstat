#
#	Kmulti.inhom.S		
#
#	$Revision: 1.7 $	$Date: 2006/06/14 14:44:57 $
#
#
# ------------------------------------------------------------------------

"Kcross.inhom" <- 
function(X, i=1, j=2, lambdaI, lambdaJ, ..., lambdaIJ=NULL)
{
	verifyclass(X, "ppp")
	if(!is.marked(X))
		stop("point pattern has no marks (no component 'marks')")
	
	I <- (X$marks == i)
	J <- (X$marks == j)
        Iname <- paste("points with mark i =", i)
        Jname <- paste("points with mark j =", j)
	result <- Kmulti.inhom(X, I, J, lambdaI, lambdaJ, ...,
                               lambdaIJ=lambdaIJ, Iname=Iname, Jname=Jname)
        attr(result, "ylab") <- substitute(Kcross.inhom(r), NULL)
        result
}

"Kdot.inhom" <- 
function(X, i=1, lambdaI, lambdadot, ..., lambdaIdot=NULL)
{
	verifyclass(X, "ppp")
	if(!is.marked(X))
		stop("point pattern has no marks (no component 'marks')")
	
	I <- (X$marks == i)
	J <- rep(TRUE, X$n)  # i.e. all points
        Iname <- paste("points with mark i =", i)
        Jname <- paste("points")
	
	result <- Kmulti.inhom(X, I, J, lambdaI, lambdadot, ...,
                         lambdaIJ=lambdaIdot,
                         Iname=Iname, Jname=Jname)
        attr(result, "ylab") <- substitute(Kdot.inhom(r), NULL)
        result
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
        correction.given <- !missing(correction) && (correction != NULL)
        correction.name <- c("border", "bord.modif", "isotropic", "Ripley", "translate")
        correction.list <- c("border", "bord.modif", "isotropic", "isotropic", "translate")
        correction.id <- pmatch(correction, correction.name)
        if(any(nbg <- is.na(correction.id)))
          stop(paste("unrecognised correction",
                     if(sum(nbg) > 1) "s",
                     ": ",
                     paste(correction[nbg], collapse=", "),
                     sep=""))
        correction <- correction.list[correction.id]
        
        # available selection of edge corrections depends on window
        if(W$type == "mask") {
           iso <- (correction == "isotropic") | (correction == "Ripley")
           if(all(iso))
             stop("Isotropic correction not implemented for binary masks")
           if(any(iso)) {
             if(correction.given)
               warning("Isotropic correction not implemented for binary masks")
             correction <- correction[!iso]
           }
        }

        # validate I, J 
	if(!is.logical(I) || !is.logical(J))
		stop("I and J must be logical vectors")
	if(length(I) != npoints || length(J) != npoints)
	     stop("The length of I and J must equal \
 the number of points in the pattern")
	
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
          lambdaI <- lambdaI[X[I]]
          if(any(is.na(lambdaI)))
            stop(paste("Pixel value of", sQuote("lambdaI"),
                       "was NA at some points of X"))
        } else if(is.vector(lambdaI) && is.numeric(lambdaI)) {
          # validate intensity vector
          if(length(lambdaI) != nI)
            stop(paste("The length of", sQuote("lambdaI"),
                       "should equal the number of", Iname))
        } else 
        stop(paste(sQuote("lambdaI"), "should be a vector or an image"))

        if(is.im(lambdaJ)) {
          # look up intensity values
          lambdaJ <- lambdaJ[X[J]]
          if(any(is.na(lambdaJ)))
            stop(paste("Pixel value of", sQuote("lambdaJ"),
                       "was NA at some points of X"))
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
        desc <- c("distance argument r", "theoretical Poisson K(r)")
        K <- fv(K, "r", substitute(Kmulti.inhom(r), NULL),
                "theo", , alim, c("r","Kpois(r)"), desc)

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
            K <- bind.fv(K, data.frame(border=Kb), "Kbord(r)",
                         "border-corrected estimate of Kmulti.inhom(r)",
                         "border")
          }
          if(any(correction == "bord.modif")) {
            Kbm <- RS$numerator/eroded.areas(W, r)
            K <- bind.fv(K, data.frame(bord.modif=Kbm), "Kbord*(r)",
                         "modified border-corrected estimate of Kmulti.inhom(r)",
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
            K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                         "translation-corrected estimate of Kmulti.inhom(r)",
                         "trans")
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(XI[icloseI], matrix(dclose, ncol=1))
            allweight <- edgewt * weight
            wh <- whist(dclose, breaks$val, allweight)
            Kiso <- cumsum(wh)/area
            rmax <- diameter(W)/2
            Kiso[r >= rmax] <- NA
            K <- bind.fv(K, data.frame(iso=Kiso), "Kiso(r)",
                         "Ripley isotropic correction estimate of K(r)",
                         "iso")
        }
        # default is to display them all
        attr(K, "fmla") <- . ~ r
        return(K)
          
}
