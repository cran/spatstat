#
#    plot.ppm.S
#
#    $Revision: 1.11 $    $Date: 2002/05/13 12:41:10 $
#
#    plot.ppm()
#         Plot a point process model fitted by mpl().
#        
#
#
plot.ppm <- function(x, nx = 40, ny = 40, 
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
        stationary <- is.stationary.ppm(model)
        poisson <- is.poisson.ppm(model)
        markeddata <- is.marked(model$Q$data)
        funnywindow <- (model$Q$data$window$type != "rectangle")

        data <- model$Q$data
        
        if(markeddata) {
          mrks <- data$marks
          if(!is.factor(mrks))
            stop("Marks are not a factor -- I don't know how to plot this")
          mrkvals <- levels(mrks)
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
          
	types <- c("trend","cif")[c(trend,cif)]
	for(ttt in types) {
          # compute the predictions
          xyz <-  predict.ppm(model, nx = nx, ny = ny, type=ttt)
          # 
          if(funnywindow) {
            # determine which prediction points are inside the window
            xy <- expand.grid(x=xyz$x, y=xyz$y)
            inside <- inside.owin(xy$x, xy$y, model$Q$data$window)
          }
          if(!markeddata) {
            # can use xyz directly
            #
            # first set predictions outside window to NA
            if(funnywindow)
              xyz$z[!inside] <- NA
            #
            for(style in how) {
              switch(style, 
                persp = persp(xyz, ...),
                image = {
                  plot(data$window, main="")
                  image(xyz, add=TRUE, ...)
                  if(superimpose) plot(data, add=TRUE)
                },
                contour = {
                  plot(data$window, main="")
                  contour(xyz, add=TRUE, ...)
                  if(superimpose) plot(data, add=TRUE)
                },
                {stop(paste("Unrecognised plot style", style))}
              )
            }
          } else {
            # marked points - slice the prediction at each level
            slice <- xyz
            for(level in mrkvals) {
              slice$z <- xyz$z[,,level]
              # set predictions outside window to NA
              if(funnywindow)
                slice$z[!inside] <- NA
              #
              main <- paste("mark =",level)
              for(style in how) {
                switch(style, 
                  persp = {
                    persp(slice, ...)
                    title(main=main)
                  },
                  image = {
                    plot(data$window, main=main)
                    image(slice, add=TRUE, ...)
                    if(superimpose) plot(data[data$marks == level], add=TRUE)
                  },
                  contour = {
                    plot(data$window, main=main)
                    contour(slice, add=TRUE, ...)
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
        if(pause)
          par(oldpar)
        return(invisible())
}      

