\name{smooth.fv}
\alias{smooth.fv}
\title{
  Apply Smoothing to Function Values
}
\description{
  Applies spline smoothing to the values
  in selected columns of a function value table.
}
\usage{
smooth.fv(x, which = ".", ...)
}
\arguments{
  \item{x}{
    Values to be smoothed.
    A function value table (object of class \code{"fv"},
    see \code{\link{fv.object}}).
  }
  \item{which}{
    Character vector identifying which columns of the table
    should be smoothed. Either a vector containing names
    of columns, or the string \code{"."} representing all the
    columns that are plotted by default.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{smooth.spline}} to
    control the smoothing.
  }
}
\details{
  \code{\link{smooth.spline}} is applied 
  to each of the selected columns in turn
  (using the function argument as the \eqn{x} coordinate
  and the selected column as the \eqn{y} coordinate).
  The original function values are then replaced by the corresponding
  smooth interpolated function values.
}
\value{
  Another function value table (object of class \code{"fv"})
  of the same format.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{with.fv}},
  \code{\link{fv.object}},
  \code{\link{smooth.spline}}
}
\examples{
   data(cells)
   G <- Gest(cells)
   plot(G)
   plot(smooth.fv(G, df=9), add=TRUE)
}
\keyword{spatial}
\keyword{nonparametric}

