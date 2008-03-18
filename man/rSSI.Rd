\name{rSSI}
\alias{rSSI}
\title{Simulate Simple Sequential Inhibition}
\description{
  Generate a random point pattern, a realisation of the
  Simple Sequential Inhibition (SSI) process.
}
\usage{
 rSSI(r, n, win = owin(c(0,1),c(0,1)), giveup = 1000)
}
\arguments{
  \item{r}{
    Inhibition distance.
  }
  \item{n}{
    Number of points to generate.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{giveup}{
    Number of rejected proposals after which the algorithm should terminate.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  This algorithm generates a realisation of the Simple Sequential
  Inhibition point process 
  inside the window \code{win}.

  Starting with an empty window, the algorithm adds points
  one-by-one. Each new point is generated uniformly in the window
  and independently of preceding points. If the new point lies
  closer than \code{r} units from an existing point, then it is rejected
  and another random point is generated.

  The algorithm terminates either when the desired number \code{n}
  of points is reached, or when the current point configuration
  has not changed for \code{giveup} iterations.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMaternI}},
\code{\link{rMaternII}}
}
\examples{
 pp <- rSSI(0.05, 200)
 \dontrun{plot(pp)}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

