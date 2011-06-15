\name{linearK}
\alias{linearK}
\title{
  Network K Function
}
\description{
  Computes an estimate of the Okabe-Yamada network \eqn{K} function
  for a point pattern on a linear network.
}
\usage{
linearK(X, r=NULL, ...)
}
\arguments{
  \item{X}{
    Point pattern on linear network (object of class \code{"lpp"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}.
    There is a sensible default.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This command computes the
  network \eqn{K} function, as defined by Okabe and Yamada (2001),
  for a point pattern on a linear network.
}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\references{
  Okabe, A. and Yamada, I. (2001) The K-function method on a network and
  its computational implementation. \emph{Geographical Analysis}
  \bold{33}, 271-290.
}

\seealso{
  \code{\link{compileK}},
  \code{\link{lpp}}
}
\examples{
  data(simplenet)
  X <- rpoislpp(5, simplenet)
  linearK(X)
}
\keyword{spatial}
\keyword{nonparametric}
