\name{is.marked.ppm}
\alias{is.marked.ppm}
\title{Test Whether A Point Process Model is Marked}
\description{
  Tests whether a fitted point process model involves ``marks''
  attached to the points.
}
\usage{
  is.marked.ppm(X, \dots) 
}
\arguments{
  \item{X}{
    Fitted point process model (object of class \code{"ppm"})
    usually obtained from \code{\link{mpl}}.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Logical value, equal to \code{TRUE} if
  \code{X} is a model that was fitted to a marked point pattern dataset.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link{longleaf}} dataset contains the locations
  of trees, each tree being marked by its diameter;
  the \code{\link{amacrine}} dataset gives the locations of cells
  of two types (on/off) and the type of cell may be regarded as a mark attached
  to the location of the cell.

  The argument \code{X} is a fitted point process model
  (an object of class \code{"ppm"}) typically obtained
  by fitting a model to point pattern data using \code{\link{mpl}}.

  This function returns \code{TRUE} if the \emph{original data}
  (to which the model \code{X} was fitted) were a marked point pattern.

  Note that this is not the same as testing whether the
  model involves terms that depend on the marks (i.e. whether the
  fitted model ignores the marks in the data).
  Currently we have not implemented a test for this.

  If this function returns \code{TRUE}, the implications are
  (for example) that
  any simulation of this model will require simulation of random marks
  as well as random point locations.
}
\seealso{
  \code{\link{is.marked}},
  \code{\link{is.marked.ppm}}
}
\examples{
  require(spatstat)
  
  data(lansing)
  # Multitype point pattern --- trees marked by species

  fit1 <- mpl(lansing, ~ marks, Poisson())
  is.marked(fit1)
  # TRUE

  fit2 <- mpl(lansing, ~ 1, Poisson())
  is.marked(fit2)
  # TRUE

  data(cells)
  # Unmarked point pattern

  fit3 <- mpl(cells, ~ 1, Poisson())
  is.marked(fit3)
  # FALSE

}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
