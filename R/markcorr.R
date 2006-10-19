#
#
#     markcorr.R
#
#     $Revision: 1.17 $ $Date: 2006/10/18 06:06:54 $
#
#    Estimate the mark correlation function
#
#
# ------------------------------------------------------------------------

"markcorr"<-
function(X, f = function(m1,m2) { m1 * m2}, r=NULL, 
         correction=c("isotropic", "Ripley", "translate"),
         method="density", ...)
{
	verifyclass(X, "ppp")
        if(!is.function(f))
          stop("Second argument f must be a function")

	npoints <- X$n
        W <- X$window

        # r values 
        rmaxdefault <- rmax.rule("K", W, npoints/area.owin(W))
        breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
        r <- breaks$r
        rmax <- breaks$max
        
        if(length(method) > 1)
          stop("Select only one method, please")
        if(method=="density" && !breaks$even)
          stop(paste("Evenly spaced r values are required if method=",
                     sQuote("density"), sep=""))
        
        # available selection of edge corrections depends on window
        if(W$type == "mask") {
          iso <- (correction == "isotropic") | (correction == "Ripley")
          if(any(iso)) {
            if(!missing(correction))
              warning("Isotropic correction not implemented for binary masks")
            correction <- correction[!iso]
          }
        }
         
        # this will be the output data frame
        result <- data.frame(r=r, theo= rep(1,length(r)))
        desc <- c("distance argument r", "theoretical Poisson m(r)=1")
        alim <- c(0, min(rmax, rmaxdefault))
        result <- fv(result,
                     "r", substitute(m(r), NULL), "theo", , alim, c("r","mpois(r)"), desc)


        # find close pairs of points
        close <- closepairs(X, rmax)
        dIJ <- close$d
        I   <- close$i
        J   <- close$j
        XI <- ppp(close$xi, close$yi, window=W)

        # apply f to each combination of marks
        #
        marx <- marks(X, dfok=FALSE)
        ff <- f(marx[I], marx[J])

        Ef <- mean(ff)
        
        if(any(correction == "translate")) {
          # translation correction
          XJ <- ppp(close$xj, close$yj, window=W)
          edgewt <- edge.Trans(XI, XJ, paired=TRUE)
          # get smoothed estimate of mark covariance
          Mtrans <- mkcor(dIJ, ff, edgewt, Ef, r, method, ...)
          result <- bind.fv(result,
                            data.frame(trans=Mtrans), "mtrans(r)",
                            "translation-corrected estimate of m(r)",
                            "trans")
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
          edgewt <- edge.Ripley(XI, matrix(dIJ, ncol=1))
          # get smoothed estimate of mark covariance
          Miso <- mkcor(dIJ, ff, edgewt, Ef, r, method, ...)
          result <- bind.fv(result,
                            data.frame(iso=Miso), "miso(r)",
                            "Ripley isotropic correction estimate of m(r)",
                            "iso")
        }
        # which corrections have been computed?
        nama2 <- names(result)
        corrxns <- rev(nama2[nama2 != "r"])

        # default is to display them all
        attr(result, "fmla") <- deparse(as.formula(paste(
                       "cbind(",
                        paste(corrxns, collapse=","),
                        ") ~ r")))
        units(result) <- units(X)
        return(result)
}
	
mkcor <- function(d, ff, wt, Ef, rvals, method="smrep", ..., nwtsteps=500) {
  d <- as.vector(d)
  ff <- as.vector(ff)
  wt <- as.vector(wt)
  switch(method,
         density={
           fw <- ff * wt
           sum.fw <- sum(fw)
           sum.wt <- sum(wt)
           # smooth estimate of kappa_f
           est <- density(d, weights=fw/sum.fw,
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           numerator <- est * sum.fw
           # smooth estimate of kappa_1
           est0 <- density(d, weights=wt/sum.wt, 
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           denominator <- est0 * Ef * sum.wt
           result <- numerator/denominator
         },
         sm={
           # This is slow!
           require(sm)
           # smooth estimate of kappa_f
           fw <- ff * wt
           est <- sm.density(d, weights=fw,
                             eval.points=rvals,
                             display="none", nbins=0, ...)$estimate
           numerator <- est * sum(fw)/sum(est)
           # smooth estimate of kappa_1
           est0 <- sm.density(d, weights=wt,
                              eval.points=rvals,
                              display="none", nbins=0, ...)$estimate
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
           result <- numerator/denominator
         },
         smrep={
           require(sm)

           hstuff <- resolve.defaults(list(...), list(hmult=1, h.weights=NA))
           if(hstuff$hmult == 1 && all(is.na(hstuff$h.weights)))
             warning("default smoothing parameter may be inappropriate")
           
           # use replication to effect the weights (it's faster)
           nw <- round(nwtsteps * wt/max(wt))
           drep.w <- rep(d, nw)
           fw <- ff * wt
           nfw <- round(nwtsteps * fw/max(fw))
           drep.fw <- rep(d, nfw)

           # smooth estimate of kappa_f
           est <- sm.density(drep.fw,
                             eval.points=rvals,
                             display="none", ...)$estimate
           numerator <- est * sum(fw)/sum(est)
           # smooth estimate of kappa_1
           est0 <- sm.density(drep.w,
                              eval.points=rvals,
                              display="none", ...)$estimate
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
           result <- numerator/denominator
         },
         loess = {
           if(!exists("R.Version") ||
              ((virg <- R.Version())$major == "1" && as.numeric(virg$minor) < 9))
             require(modreg)
           # set up data frame
           df <- data.frame(d=d, ff=ff, wt=wt)
           # fit curve to numerator using loess
           fitobj <- loess(ff ~ d, data=df, weights=wt, ...)
           # evaluate fitted curve at desired r values
           Eff <- predict(fitobj, newdata=data.frame(d=rvals))
           # normalise:
           # denominator is the sample mean of all ff[i,j],
           # an estimate of E(ff(M1,M2)) for M1,M2 independent marks
           result <- Eff/Ef
         },
         )
  return(result)
}

