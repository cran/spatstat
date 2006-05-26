\name{rmh}
\alias{rmh}
\title{Simulate point patterns using the Metropolis-Hastings algorithm.}
\description{
  Generic function for running the Metropolis-Hastings algorithm
  to produce simulated realisations of a point process model.
}

\usage{rmh(model, \dots)}

\arguments{
  \item{model}{The point process model to be simulated.
  }
  \item{\dots}{Further arguments controlling the simulation.
  }
}

\details{
  The Metropolis-Hastings algorithm can be used to
  generate simulated realisations from a wide range of
  spatial point processes. For caveats, see below.
  
  The function \code{rmh} is generic; it has methods
  \code{\link{rmh.ppm}} (for objects of class \code{"ppm"})
  and  \code{\link{rmh.default}} (the default).
  The actual implementation of the Metropolis-Hastings algorithm is
  contained in \code{\link{rmh.default}}.
  For details of its use, see 
  \code{\link{rmh.ppm}} or \code{\link{rmh.default}}.

  [If the model is a Poisson process, then Metropolis-Hastings
  is not used; the Poisson model is generated directly
  using \code{\link{rpoispp}} or \code{\link{rmpoispp}}.]

  In brief, the Metropolis-Hastings algorithm is a Markov Chain,
  whose states are spatial point patterns, and whose limiting
  distribution is the desired point process. After
  running the algorithm for a very large number of iterations,
  we may regard the state of the algorithm as a realisation
  from the desired point process.

  However, there are difficulties in deciding whether the
  algorithm has run for ``long enough''. The convergence of the
  algorithm may indeed be extremely slow. No guarantees of
  convergence are given!

  While it is fashionable to decry the Metropolis-Hastings algorithm
  for its poor convergence and other properties, it has the advantage
  of being easy to implement for a wide range of models.
}

\value{
  A point pattern, in the form of an object of class \code{"ppp"}.
  See \code{\link{rmh.default}} for details.
}

\seealso{
  \code{\link{rmh.default}}
}

\examples{
    # See examples in rmh.default and rmh.ppm
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

