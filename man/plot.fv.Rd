\name{plot.fv}
\alias{plot.fv}
\title{Plot Function Valuesn}
\description{
  Plot method for the class \code{"fv"}.
}
\usage{
plot.fv(x, fmla, subset=NULL, lty=NULL, col=NULL, xlim, ylim, xlab, ylab, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{"fv"}, containing the variables to be plotted
    or variables from which the plotting coordinates can be computed.
  }
  \item{fmla}{
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
  \item{\dots}{
    Extra arguments passed to \code{plot.default}.
  }
}
\value{
  none. Side effect is a plot.
}
\details{
  This is the \code{plot} method for the class \code{"fv"}.

  The use of the argument \code{fmla} is like \code{plot.formula}, but offers
  some extra functionality.
  
  The left and right hand sides of \code{fmla}
  are evaluated in the data frame \code{x},
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
   data(cells)
   K <- Kest(cells)
   # K is an object of class "fv"

   plot(K, iso ~ r)                # plots iso against r

   plot(K, sqrt(iso/pi) ~ r)   # plots sqrt(iso/r)  against r

   plot(K, cbind(iso,theo) ~ r)   # plots iso against r  AND theo against r

   plot(K, cbind(iso,theo) ~ r, col=c(2,3))
                                   # plots iso against r  in colour 2
                                   # and theo against r in colour 3

   plot(K, iso ~ r, subset=quote(r < 0.2))
                                   # plots iso against r for r < 10
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{Kest}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
