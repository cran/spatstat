\name{area.owin}
\alias{area.owin}
\title{Area of a Window}
\description{
  Computes the area of a window 
}
\usage{
 area.owin(w)
}
\arguments{
  \item{w}{A window, whose area will be computed.
    This should be an object of class \code{\link{owin}},
    or can be given in any format acceptable to \code{\link{as.owin}()}.
  }
}
\value{
  The numerical value of the area of the window. 
}
\details{
  If the window \code{W} is of type \code{"rectangle"} or \code{"polygonal"},
  the area of this rectangular window is computed by analytic geometry.
  If \code{W} is of type \code{"mask"} 
  the area of the discrete raster approximation of the window is
  computed by summing the binary image values and adjusting for
  pixel size.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{as.owin}}
}
\examples{
  w <- unit.square()
  area.owin(w)
       # returns 1.00000

  k <- 6
  theta <- 2 * pi * (0:(k-1))/k
  co <- cos(theta)
  si <- sin(theta)
  mas <- owin(c(-1,1), c(-1,1), poly=list(x=co, y=si))
  area.owin(mas)
      # returns approx area of k-gon
  
  mas <- as.mask(w, eps=0.01)
  X <- raster.x(mas)
  Y <- raster.y(mas)
  mas$m <- ((X - 0.5)^2 + (Y - 0.5)^2 <= 1)
  area.owin(mas)
       # returns 3.14 approx     

}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{math}
