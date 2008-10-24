\name{quadratcount}
\alias{quadratcount}
\title{Quadrat counting for a point pattern}
\description{
  Divides window into quadrats and 
  counts the numbers of points in each quadrat.
}
\usage{
  quadratcount(X, nx=5, ny=nx, \dots,
               xbreaks=NULL, ybreaks=NULL, tess=NULL)
}
\arguments{
  \item{X}{
    A point pattern
    (object of class \code{"ppp"}).
    }
    \item{nx,ny}{
      Numbers of rectangular quadrats in the \eqn{x} and \eqn{y} directions.
      Incompatible with \code{xbreaks} and \code{ybreaks}.
    }
    \item{\dots}{Ignored.}
    \item{xbreaks}{
      Numeric vector giving the \eqn{x} coordinates of the
      boundaries of the rectangular quadrats. Incompatible with \code{nx}.
    }
    \item{ybreaks}{
      Numeric vector giving the \eqn{y} coordinates of the
      boundaries of the rectangular quadrats. Incompatible with \code{ny}.
    }
    \item{tess}{
      Tessellation (object of class \code{"tess"})
      determining the quadrats. Incompatible
      with \code{nx,ny,xbreaks,ybreaks}.
    }
}
\value{
  A contingency table containing the number of points in each
  quadrat.

  The table is also an object of the special class \code{"quadratcount"}
  and there is a plot method for this class.
}
\details{
  Quadrat counting is an elementary technique for analysing spatial
  point patterns. See Diggle (2003).

  By default, the window containing the point pattern \code{X} is divided into
  an \code{nx * ny} grid of rectangular tiles or `quadrats'.
  (If the window is not a rectangle, then these tiles are intersected
  with the window.)
  The number of points of \code{X} falling in each quadrat is
  counted. These numbers are returned as a contingency table.

  If \code{xbreaks} is given, it should be a numeric vector
  giving the \eqn{x} coordinates of the quadrat boundaries.
  If it is not given, it defaults to a
  sequence of \code{nx+1} values equally spaced
  over the range of \eqn{x} coordinates in the window \code{X$window}.

  Similarly if \code{ybreaks} is given, it should be a numeric
  vector giving the \eqn{y} coordinates of the quadrat boundaries.
  It defaults to a vector of \code{ny+1} values
  equally spaced over the range of \eqn{y} coordinates in the window.
  The lengths of \code{xbreaks} and \code{ybreaks} may be different.

  Alternatively, quadrats of any shape may be used.
  The argument \code{tess} can be a tessellation (object of class
  \code{"tess"}) whose tiles will serve as the quadrats.
  
  The algorithm counts the number of points of \code{X}
  falling in each quadrat, and returns these counts as a
  contingency table. 

  The return value is a \code{table} which can be printed neatly.
  The return value is also a member of the special class
  \code{"quadratcount"}. Plotting the object will display the
  quadrats, annotated by their counts. See the examples.
  
  To perform a chi-squared test based on the quadrat counts,
  use \code{\link{quadrat.test}}.
}
\seealso{
  \code{\link{quadrat.test}},
  \code{\link{quadrats}},
  \code{\link{quadratresample}}
}
\references{
  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Academic Press, 2003.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\examples{
 X <- runifpoint(50)
 quadratcount(X)
 quadratcount(X, 4, 5)
 quadratcount(X, xbreaks=c(0, 0.3, 1), ybreaks=c(0, 0.4, 0.8, 1))
 qX <-  quadratcount(X, 4, 5)

 # plotting:
 plot(X, pch="+")
 plot(qX, add=TRUE, col="red", cex=1.5, lty=2)

 # quadrats determined by tessellation:
 B <- dirichlet(runifpoint(6))
 qX <- quadratcount(X, tess=B)
 plot(X, pch="+")
 plot(qX, add=TRUE, col="red", cex=1.5, lty=2)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
 
 
