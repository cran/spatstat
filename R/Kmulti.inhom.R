#
#	Kmulti.inhom.S		
#
#	$Revision: 1.1 $	$Date: 2005/10/14 09:40:20 $
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

        breaks <- handle.r.b.args(r, breaks, W)
        r <- breaks$r

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
        if(W$type != "rectangle") {
           iso <- (correction == "isotropic") | (correction == "Ripley")
           if(all(iso))
             stop("Isotropic correction not implemented for non-rectangular windows")
           if(any(iso)) {
             if(correction.given)
               warning("Isotropic correction not implemented for non-rectangular windows")
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

        # validate intensity vectors
        if(!is.vector(lambdaI))
          stop(paste(sQuote("lambdaI"), "should be a vector"))
        if(length(lambdaI) != nI)
          stop(paste("The length of", sQuote("lambdaI"),
                     "should equal the number of", Iname))

        if(!is.vector(lambdaJ))
          stop(paste(sQuote("lambdaJ"), "should be a vector"))
        if(length(lambdaJ) != nJ)
          stop(paste("The length of", sQuote("lambdaJ"),
                     "should equal the number of", Jname))

        # Form weight for each pair
        if(is.null(lambdaIJ))
          weight <- 1/outer(lambdaI, lambdaJ, "*")
        else {
          if(!is.matrix(lambdaIJ))
            stop("lambdaIJ should be a matrix")
          if(nrow(lambdaIJ) != nI)
            stop(paste("nrow(lambdaIJ) should equal the number of", Iname))
          if(ncol(lambdaIJ) != nJ)
            stop(paste("ncol(lambdaIJ) should equal the number of", Jname))
          weight <- 1/lambdaIJ
        }

        # Recommended range of r values
        alim <- c(0, min(diff(X$window$xrange), diff(X$window$yrange))/4)
        
        # this will be the output data frame
        # It will be given more columns later
        K <- data.frame(r=r, theo= pi * r^2)
        desc <- c("distance argument r", "theoretical Poisson K(r)")
        K <- fv(K, "r", substitute(Kmulti.inhom(r), NULL),
                "theo", , alim, c("r","Kpois(r)"), desc)

# interpoint distances		
	d <- crossdist(x[I], y[I], x[J], y[J])
# distances to boundary	
	b <- (bdist.points(X))[I]
        
# Determine which interpoint distances d[i,j] refer to the same point
# (not just which distances are zero)        
        same <- matrix(FALSE, nrow=nI, ncol=nJ)
        common <- I & J
        if(any(common)) {
          Irow <- cumsum(I)
          Jcol <- cumsum(J)
          icommon <- (1:npoints)[common]
          for(i in icommon)
            same[Irow[i], Jcol[i]] <- TRUE
        }

# Compute estimates by each of the selected edge corrections.
        
        if(any(correction == "border" | correction == "bord.modif")) {
          # border method
          # Compute distances to boundary
          b <- bdist.points(X[I])
          # Distances corresponding to identical pairs
          # are excluded from consideration
          d[same] <- Inf
          # apply reduced sample algorithm
          RS <- Kwtsum(d, b, weight, 1/lambdaI, breaks, slow=FALSE)
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
          # reset identical pairs to original values
          d[same] <- 0
        }
        if(any(correction == "translate")) {
          # translation correction
            edgewt <- edge.Trans(X[I], X[J])
            allweight <- edgewt * weight
            wh <- whist(d[!same], breaks$val, allweight[!same])
            Ktrans <- cumsum(wh)/area
            rmax <- diameter(W)/2
            Ktrans[r >= rmax] <- NA
            K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                         "translation-corrected estimate of Kmulti.inhom(r)",
                         "trans")
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(X[I], d)
            allweight <- edgewt * weight
            wh <- whist(d[!same], breaks$val, allweight[!same])
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
