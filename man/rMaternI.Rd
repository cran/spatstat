\name{rMaternI}
\alias{rMaternI}
\title{Simulate Matern Model I}
\description{
  Generate a random point pattern using the
  Mat\'ern Model I inhibition process.
}
\usage{
 rMaternI(kappa, r, win = owin(c(0,1),c(0,1)))
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of proposal points.
    A single positive number.
  }
  \item{r}{
    Inhibition distance.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  This algorithm generates a realisation of Mat\'ern's Model I
  inhibition process 
  inside the window \code{win}.

  The process is constructed by first
  generating a uniform Poisson point process of ``proposal'' points 
  with intensity \code{kappa} inside the window.
  A proposal point is deleted if it lies within \code{r} units' distance
  of another proposal point. Otherwise it is retained.
  The retained points constitute Mat\'ern's Model I.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}}
}
\examples{
 pp <- rMaternI(20, 0.05)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{datagen}

