\name{rpoispp}
\alias{rpoispp}
\title{Generate Poisson Point Pattern}
\description{
  Generate a random point pattern using the
  (homogeneous or inhomogeneous) Poisson process.
}
\usage{
 rpoispp(lambda, max, win, \dots)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    Either a single positive number, or a \code{function(x,y, \dots)}.
  }
  \item{max}{
    An upper bound for the value of \code{lambda(x,y)},
    if \code{lambda} is a function.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{\dots}{
    Arguments passed to \code{lambda} if it is a function.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  If \code{lambda} is a single number,
  then this algorithm generates a realisation
  of the uniform Poisson process inside the window \code{win} with 
  intensity \code{lambda} (points per unit area).
 
  If \code{lambda} is a function, then this algorithm generates a realisation
  of the inhomogeneous Poisson process with intensity function
  \code{lambda(x,y,\dots)} at spatial location \code{(x,y)}.
  The function \code{lambda} must work correctly with vectors \code{x}
  and \code{y}.
  The value \code{max} must be given and must be an upper bound on the 
  values of \code{lambda(x,y,\dots)} for all locations \code{(x, y)}
  inside the window \code{win}. 
 
  To generate an inhomogeneous Poisson process
  the algorithm uses ``rejection filtering'': it first generates a uniform
  Poisson process of intensity \code{max},
  then thins it by randomly deleting or retaining each point independently,
  with retention probability
  \eqn{p(x,y) = \lambda(x,y)/\mbox{max}}{p(x,y) = lambda(x,y)/max}.
}
\seealso{
\code{\link{ppp.object}},
\code{\link{owin.object}}
}
\examples{
 # uniform Poisson process with intensity 100 in the unit square
 pp <- rpoispp(100)
 
 # uniform Poisson process with intensity 1 in a 10 x 10 square
 pp <- rpoispp(1, win=owin(c(0,10),c(0,10)))
 # plots should look similar !
 
 # inhomogeneous Poisson process in unit square
 # with intensity lambda(x,y) = 100 * exp(-3*x)
 # Intensity is bounded by 100
 pp <- rpoispp(function(x,y) {100 * exp(-3*x)}, 100)

 # How to tune the coefficient of x
 lamb <- function(x,y,a) { 100 * exp( - a * x)}
 pp <- rpoispp(lamb, 100, a=3)

}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
