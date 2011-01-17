#
#	$Revision: 1.22 $	$Date: 2011/01/17 02:46:06 $
#
#    ppm()
#          Fit a point process model to a two-dimensional point pattern
#
#

"ppm" <- 
function(Q,
         trend = ~1,
	 interaction = Poisson(),
         ..., 
         covariates = NULL,
         covfunargs = list(),
	 correction="border",
	 rbord = reach(interaction),
         use.gam=FALSE,
         method = "mpl",
         forcefit=FALSE,
         gcontrol=list(),
         nsim=100,
         nrmh=1e5,
         start=NULL,
         control=list(nrep=nrmh),
         verb=TRUE
) {
  if(!(method %in% c("mpl", "ho")))
    stop(paste("Unrecognised fitting method", sQuote(method)))
  cl <- match.call()
  callstring <- paste(deparse(sys.call()), collapse="")

  if(is.null(interaction))
    interaction <- Poisson()

  if(is.ppp(Q) && is.marked(Q) && !is.multitype(Q)) 
    stop(paste("ppm is not yet implemented for marked point patterns,",
               "other than multitype patterns."))
  
  # validate choice of edge correction
  correction <- pickoption("correction", correction,
                           c(border="border",
                             periodic="periodic",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             translation="translate",
                             none="none"))
  
  # validate rbord for border correction
  if(correction == "border") {
    rbord.given <- !missing(rbord) && !is.null(rbord)
    if(is.null(rbord))
      rbord <- reach(interaction)
    infin <- is.infinite(rbord)
    too.large <- infin || (eroded.areas(as.owin(Q), rbord) == 0)
    if(too.large) {
      whinge <-
        paste(if(rbord.given) "rbord" else "the reach of this interaction",
              if(infin) "is infinite or unknown;"
              else "is too large for this window;",
              "please specify",
              if(rbord.given) "a smaller value of",
              "rbord, or use a different edge correction")
      stop(whinge)
    }
  }
  
  # go
  fitMPL <- mpl.engine(Q=Q, trend=trend,
                       interaction=interaction,
                       covariates=covariates,
                       covfunargs=covfunargs,
                       correction=correction,
                       rbord=rbord,
                       use.gam=use.gam,
                       forcefit=forcefit,
                       gcontrol=gcontrol,
                       callstring=callstring,
                       preponly=FALSE,
                       ...)
  
  fitMPL$call <- cl
  fitMPL$callstring <- callstring
  fitMPL$callframe <- parent.frame()

  if(method == "mpl" || is.poisson.ppm(fitMPL))
    return(fitMPL)

  fitHO <- ho.engine(fitMPL, nsim=nsim, nrmh=nrmh, start=start,
                     control=control, verb=verb)

  if(!is.null(fitHO))
    return(fitHO)
  else
    return(fitMPL)

}

