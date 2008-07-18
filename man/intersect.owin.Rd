\name{intersect.owin}
\alias{intersect.owin}
\alias{union.owin}
\title{Intersection or Union of Two Windows}
\description{
  Yields the intersection or union of two windows.
}
\usage{
   intersect.owin(A, B, \dots)
   union.owin(A,B, \dots)
}
\arguments{
  \item{A}{A window object (see Details).}
  \item{B}{A window object.}
  \item{\dots}{Optional arguments passed to \code{\link{as.mask}}
    to control the discretisation, if required.
    }
}
\value{
  A window (object of class \code{"owin"}).
}
\details{
  The function \code{intersect.owin} computes the intersection between the
  two windows \code{A} and \code{B}, while
  \code{union.owin} computes their union.

  The arguments \code{A} and \code{B} must be window objects
  (either objects of class \code{"owin"}, or data that can be
  coerced to this class by \code{\link{as.owin}}).

  If the intersection or union is empty, an error occurs.
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{is.subset.owin}},
  \code{\link{overlap.owin}},
  \code{\link{bounding.box}},
  \code{\link{owin.object}}
}
\examples{
# rectangles
   u <- unit.square()
   v <- owin(c(0.5,3.5), c(0.4,2.5))
# polygon
   data(letterR)
# mask
   m <- as.mask(letterR)

# two rectangles
   intersect.owin(u, v) 
   union.owin(u,v)

# polygon and rectangle
   intersect.owin(letterR, v)
   union.owin(letterR,v)

# mask and rectangle
   intersect.owin(m, v)
   union.owin(m,v)

# mask and polygon
   p <- rotate(v, 0.2)
   intersect.owin(m, p)
   union.owin(m,p)

# two polygons
   A <- letterR
   B <- rotate(letterR, 0.2)
   plot(bounding.box(A,B), main="intersection")
   w <- intersect.owin(A, B)
   plot(w, add=TRUE)
   plot(A, add=TRUE)
   plot(B, add=TRUE)

   plot(bounding.box(A,B), main="union")
   w <- union.owin(A,B)
   plot(w, add=TRUE)   
   plot(A, add=TRUE)
   plot(B, add=TRUE)
}
\keyword{spatial}
\keyword{math}
