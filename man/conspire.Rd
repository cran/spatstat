\name{conspire}
\alias{conspire}
\title{Plot Together According to a Formula}
\description{
  Given a data frame and a formula involving variables in the data frame,
  plot the left hand side of the formula versus the right hand side. 
}
\usage{
 conspire(indata, formula, subset=NULL, lty=NULL, col=NULL, xlim, ylim,
xlab, ylab)
}
\arguments{
  \item{indata}{
    A data frame, containing the variables to be plotted
    or variables from which the plotting coordinates can be computed.
  }
  \item{formula}{
    an S language formula determining which variables or expressions
    are plotted. 
  }
  \item{subset}{
    (optional) subset of rows of the data frame that will be plotted.
  }
  \item{lty}{
    (optional) numeric vector of values of the graphical parameter
    \code{lty} controlling the line style of each plot.
  }
  \item{col}{
    (optional) numeric vector of values of the graphical parameter
    \code{col} controlling the colour of each plot.
  }
  \item{xlim}{
    (optional) range of x axis
  }
  \item{ylim}{
    (optional) range of y axis
  }
  \item{xlab}{
    (optional) label for x axis
  }
  \item{ylab}{
    (optional) label for y axis
  }
}
\value{
  none. Side effect is a plot.
}
\details{
  This function is like \code{plot.formula} but offers
  some extra functionality.
  
  The left and right hand sides of the formula argument
  are evaluated in the data frame \code{indata},
  and the results are plotted against each other
  (the left side on the \eqn{y} axis 
  against the right side on the \eqn{x} axis). 
  Both left and right sides may be variables in the data frame
  or expressions in these variables.

  Multiple curves may be specified by a single formula
  of the form 
  \code{cbind(y1,y2,\ldots,yn) ~ x}, where \code{x,y1,y2,\ldots,yn} are
  expressions involving the variables in the data frame.
  Each of the variables \code{y1,y2,\ldots,yn} in turn will be plotted
  against \code{x}. 
  See the examples.
}
\examples{
   \dontrun{
   library(spatstat)

   df <- data.frame(x=1:20, y=sqrt(1:20), z = sqrt(1:20) + runif(20))

   conspire(df, y ~ x)             # plots y against x

   conspire(df, log(y) ~ log(x))   # plots log(y) against log(x)

   conspire(df, z - sqrt(x) ~ x)   # plots z - sqrt(x) against x


   conspire(df, cbind(y,z) ~ x)   # plots y against x  AND z against x

   conspire(df, cbind(z,sqrt(x)) ~ x)
                                  # plots z against x  AND sqrt(x) against x

   conspire(df, cbind(z - sqrt(x),0) ~ x)
                                   # plots z - sqrt(x) against x
                                   # with dashed horizontal line at height 0

   conspire(df, cbind(y,z) ~ x, col=c(2,3))
                                   # plots y against x  in colour 2
                                   # and z against x in colour 3

   conspire(df, y ~ x, subset=quote(x < 10))
                                   # plots y against x for x < 10

  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
