#
#	plot.ppp.S
#
#	$Revision: 1.7 $	$Date: 2002/01/18 06:39:34 $
#
#
#--------------------------------------------------------------------------

plot.ppp <- function(x, main, ..., chars, use.marks=T, add=F)
{
#
# Function plot.ppp.
# A plot() method for the class 'ppp'
#
	if(missing(main))
		main <- deparse(substitute(x))

        x <- as.ppp(x)

        if(!add)
          plot.owin(x$window, ..., main=main)

        if(x$n == 0)
          return(invisible())
        
	if(is.null(x$marks) || !use.marks) {
		points(x$x, x$y, ...)
                return(invisible())
        } else {
	        um <- if(is.factor(x$marks))
                         levels(x$marks)
                      else
                         sort(unique(x$marks))

                if(missing(chars))
                  chars <- seq(um)
                
                for(i in seq(um)) {
                  relevant <- (x$marks == um[i])
                  if(any(relevant)) 
                    points(x$x[relevant], x$y[relevant], pch = chars[i])
                }
                names(chars) <- um
                return(chars)
	}
}
