#
#    plot.ppm.S
#
#    $Revision: 1.15 $    $Date: 2004/01/08 16:21:04 $
#
#    plot.ppm()
#         Plot a point process model fitted by mpl().
#        
#
#
plot.ppm <- function(x, ngrid = c(40,40),
		     superimpose = TRUE,
                     trend = TRUE, cif = TRUE, pause = TRUE,
                     how=c("persp","image", "contour"),  ...)
{
  
        model <- x
        
#       Plot a point process model fitted by mpl().
#
        verifyclass(model, "ppm")
#
#       find out what kind of model it is
#
        mod <- summary(model)
        stationary <- mod$stationary
        poisson    <- mod$poisson
        marked     <- mod$marked
        multitype  <- mod$multitype
        funnywindow <- (mod$entries$data$window$type != "rectangle")

        data <- mod$entries$data
        
        if(marked) {
          if(!multitype)
            stop("Not implemented for general marked point processes")
          else
            mrkvals <- levels(data$marks)
        }
        
#
# Interpret options
#
#        
#                        plotting style
#        
        howmat <- outer(how, c("persp","image", "contour"), "==")
        howmatch <- apply(howmat, 1, any)
        if(any(!howmatch))
          stop(paste("unrecognised option", how[!howmatch]))

#                       whether to plot trend, cif
        
        if(!missing(trend) && !missing(cif) && !trend && !cif) {
          cat("Nothing plotted - both \'trend\' and \'cif\' are FALSE\n")
          return(invisible(NULL))
        }
#                       suppress uninteresting plots
#                       unless explicitly instructed otherwise
        if(missing(trend))
          trend <- !stationary
	if(missing(cif))
          cif <- !poisson
        
        if(!trend && !cif) {
          cat("Nothing plotted -- all plots selected are flat surfaces.\n")
          return(invisible())
        }

#
#        
#
# Do the plotting --- first the trend, then the full cif
# (if both are requested).
#

        if(pause)
          oldpar <- par(ask=TRUE)

        on.exit(if(pause) par(oldpar))

#
#        
	types <- c("trend","cif")[c(trend,cif)]
	for(ttt in types) {
          # compute the predictions
          pred <-  predict.ppm(model, ngrid=ngrid, type=ttt)
          if(!marked) {
            # can use xyz directly
            for(style in how) {
              switch(style, 
                persp = {
                  arg1 <- list(pred, xlab="x")
                  arg2 <- spatstat.options("par.persp")[[1]]
                  arg3 <- list(...)
                  argue <- append(arg1, append(arg2, arg3))
                  do.call("persp", argue)
                },
                image = {
                  image(pred, ...)
                  if(superimpose) plot(data, add=TRUE)
                },
                contour = {
                  contour(pred, ...)
                  if(superimpose) plot(data, add=TRUE)
                },
                {stop(paste("Unrecognised plot style", style))}
              )
            }
          } else {
            # multitype process - list of images
            for(i in seq(mrkvals)) {
              level <- mrkvals[i]
              main <- paste("mark =",level)
              for(style in how) {
                switch(style, 
                  persp = {
                    arg1 <- list(pred[[i]], xlab="x")
                    arg2 <- spatstat.options("par.persp")[[1]]
                    arg3 <- list(...)
                    argue <- append(arg1, append(arg2, arg3))
                    do.call("persp", argue)
                    title(main=main)
                  },
                  image = {
                    image(pred[[i]], ...)
                    if(superimpose) plot(data[data$marks == level], add=TRUE)
                  },
                  contour = {
                    contour(pred[[i]], ...)
                    if(superimpose) plot(data[data$marks == level], add=TRUE)
                  },
                  {
                    stop(paste("Unrecognised plot style", style))
                  }
                )
              }
            }
          }
        }

        return(invisible())
}      

