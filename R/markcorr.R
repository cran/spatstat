#
#
#     markcorr.R
#
#     $Revision: 1.3 $ $Date: 2004/01/13 11:33:46 $
#
#    Estimate the mark correlation function
#
#
# ------------------------------------------------------------------------

"markcorr"<-
function(X, f = function(m1,m2) { m1 * m2}, r=NULL, slow=FALSE,
         correction=c("border", "isotropic", "Ripley", "translate"),
         method="density", ...)
{
	verifyclass(X, "ppp")
        if(!is.function(f))
          stop("Second argument f must be a function")

	npoints <- X$n
        W <- X$window

        breaks <- handle.r.b.args(r, NULL, W)
        r <- breaks$r

        if(length(method) > 1)
          stop("Select only one method, please")
        if(method=="density" && !breaks$even)
          stop("Evenly spaced r values are required if method=\"density\"")
        
        # available selection of edge corrections depends on window
        if(W$type != "rectangle") {
          iso <- (correction == "isotropic") | (correction == "Ripley")
          if(any(iso)) {
            if(!missing(correction))
              warning("Isotropic correction not implemented for non-rectangular windows")
            correction <- correction[!iso]
          }
        }
         
        # this will be the output data frame
        result <- data.frame(r=r, theo= rep(1,length(r)))
        desc <- c("distance argument r", "theoretical Poisson m(r)=1")
        alim <- c(0, min(diff(X$window$xrange),diff(X$window$yrange))/4)
        result <- fv(result,
                     "r", "m(r)", "theo", , alim, c("r","mpois(r)"), desc)


        # pairwise distances
	d <- pairdist(X$x, X$y)
        offdiag <- (row(d) != col(d))

        # apply f to each combination of marks
        # mm[i,j] = mark[i]
        mm <- matrix(X$marks,nrow=npoints,ncol=npoints)
        mr <- as.vector(mm)
        mc <- as.vector(t(mm))
        ff <- f(mr, mc)
        ff <- matrix(ff, nrow=npoints, ncol=npoints)
        # so ff[i,j] = f(mark[i], mark[j])

        Ef <- mean(ff)
        
        if(any(correction == "border")) {
          # border method
          # Compute distances to boundary
          b <- bdist.points(X)
          bb <- matrix(b, nrow=npoints, ncol=npoints)
          # set edge correction weight = 1 iff b > d
          ok <- 1 * (bb > d)
          # get smoothed estimate of mcf
          Mborder <- mkcor(d[offdiag], ff[offdiag], ok[offdiag],
                                 Ef, r, method, ...)
          result <- bind.fv(result,
                            data.frame(border=Mborder), "mbord(r)",
                            "border-corrected estimate of m(r)",
                            "border")
        }
        if(any(correction == "translate")) {
          # translation correction
            edgewt <- edge.Trans(X, exact=slow)
          # get smoothed estimate of mark covariance
            Mtrans <- mkcor(d[offdiag], ff[offdiag], edgewt[offdiag],
                                  Ef, r, method, ...)
            result <- bind.fv(result,
                              data.frame(trans=Mtrans), "mtrans(r)",
                              "translation-corrected estimate of m(r)",
                              "trans")
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(X, d)
          # get smoothed estimate of mark covariance
            Miso <- mkcor(d[offdiag], ff[offdiag], edgewt[offdiag],
                                Ef, r, method, ...)
            result <- bind.fv(result,
                              data.frame(iso=Miso), "miso(r)",
                              "Ripley isotropic correction estimate of m(r)",
                              "iso")
        }
        return(result)
}
	
mkcor <- function(d, ff, wt, Ef, rvals, method="smrep", ..., nwtsteps=500) {
  d <- as.vector(d)
  ff <- as.vector(ff)
  wt <- as.vector(wt)
  switch(method,
         density={
           # use replication to effect the weights
           # (there's no weights argument to density()
           nw <- round(nwtsteps * wt/max(wt))
           drep.w <- rep(d, nw)
           fw <- ff * wt
           nfw <- round(nwtsteps * fw/max(fw))
           drep.fw <- rep(d, nfw)
           # smooth estimate of kappa_f
           est <- density(drep.fw,
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           numerator <- est * sum(fw)/sum(est)
           # smooth estimate of kappa_1
           est0 <- density(drep.w,
                          from=min(rvals), to=max(rvals), n=length(rvals),
                          ...)$y
           denominator <- est0 * (sum(wt)/ sum(est0)) * Ef
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

