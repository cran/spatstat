\name{rpoispp}
\alias{rpoispp}
\title{Generate Poisson Point Pattern}
\description{
  Generate a random point pattern using the
  (homogeneous or inhomogeneous) Poisson process.
  Includes CSR (complete spatial randomness).
}
\usage{
 rpoispp(lambda, lmax, win, \dots)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process.
    Either a single positive number, a \code{function(x,y, \dots)},
    or a pixel image.
  }
  \item{lmax}{
    An upper bound for the value of \code{lambda(x,y)},
    if \code{lambda} is a function.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
    Ignored if \code{lambda} is a pixel image.
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
  of the uniform Poisson process (also known as 
  Complete Spatial Randomness, CSR) inside the window \code{win} with 
  intensity \code{lambda} (points per unit area).
 
  If \code{lambda} is a function, then this algorithm generates a realisation
  of the inhomogeneous Poisson process with intensity function
  \code{lambda(x,y,\dots)} at spatial location \code{(x,y)}
  inside the window \code{win}.
  The function \code{lambda} must work correctly with vectors \code{x}
  and \code{y}.
  The value \code{lmax} must be given and must be an upper bound on the 
  values of \code{lambda(x,y,\dots)} for all locations \code{(x, y)}
  inside the window \code{win}. 

  If \code{lambda} is a pixel image object of class \code{"im"}
  (see \code{\link{im.object}}), this algorithm generates a realisation
  of the inhomogeneous Poisson process with intensity equal to the
  pixel values of the image. (The value of the intensity function at an
  arbitrary location is the pixel value of the nearest pixel.)
  The argument \code{win} is ignored;
  the window of the pixel image is used instead.
  
  To generate an inhomogeneous Poisson process
  the algorithm uses ``thinning'': it first generates a uniform
  Poisson process of intensity \code{lmax},
  then randomly deletes or retains each point, independently of other points,
  with retention probability
  \eqn{p(x,y) = \lambda(x,y)/\mbox{lmax}}{p(x,y) = lambda(x,y)/lmax}.
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

 # pixel image
 Z <- as.im(function(x,y){100 * sqrt(x+y)}, unit.square())
 pp <- rpoispp(Z)

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
