\name{rMaternII}
\alias{rMaternII}
\title{Simulate Matern Model II}
\description{
  Generate a random point pattern using the
  Mat\'ern Model II inhibition process.
}
\usage{
 rMaternII(kappa, r, win = owin(c(0,1),c(0,1)))
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
  This algorithm generates a realisation of Mat\'ern's Model II
  inhibition process 
  inside the window \code{win}.

  The process is constructed by first
  generating a uniform Poisson point process of ``proposal'' points 
  with intensity \code{kappa} inside the window.
  Then each proposal point is marked by an ``arrival time'', a number
  uniformly distributed in \eqn{[0,1]} independently of other variables.
  
  A proposal point is deleted if it lies within \code{r} units' distance
  of another proposal point \emph{that has an earlier arrival time}.
  Otherwise it is retained.
  The retained points constitute Mat\'ern's Model II.

  The difference between Mat\'ern's Model I and II is the italicised
  statement above. Model II has a higher intensity
  for the same parameter values.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}},
\code{\link{rMaternI}}
}
\examples{
 pp <- rMaternII(20, 0.05)
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

