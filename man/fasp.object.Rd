\name{fasp.object}
\alias{fasp.object}
\title{Function Arrays for Spatial Patterns}
\description{
  A class \code{"fasp"} to represent a ``matrix''
  of functions, amenable to plotting as a matrix of plot panels.
}
\details{
  An object of this class is a convenient way of storing
  (and later plotting, editing, etc)
  a set of functions \eqn{f_{i,j}(r)}{f[i,j](r)} of a real argument \eqn{r},
  defined for each possible pair \eqn{(i,j)} of indices
  \eqn{1 \le i,j \le n}{1 <= i,j <= n}. We may think of this
  as a matrix or array of functions \eqn{f_{i,j}}{f[i,j]}.

  Function arrays are particularly useful in the 
  analysis of a multitype point pattern (a point pattern in which 
  the points are identified as belonging to separate types).
  We may want to compute a summary function for the points
  of type \eqn{i} only, for each of the possible types \eqn{i}.
  This produces a \eqn{1 \times m}{1 * m} array of functions.
  Alternatively we may compute a summary function
  for each possible pair of types \eqn{(i,j)}.
  This produces an \eqn{m \times m}{m *  m} array of functions.

  For multitype point patterns the command \code{\link{alltypes}}
  will compute arrays of summary functions for each possible
  type or for each possible pair of types.
  For univariate (single-type) point patterns the command
  \code{\link{allstats}} will compute an array of
  different summary functions \eqn{F}, \eqn{G}, \eqn{J},
  \eqn{K} for the same dataset.
  Both \code{\link{alltypes}} and \code{\link{allstats}}
  return an object of class \code{"fasp"}.

  There are methods for \code{plot}, \code{print} and \code{"["} in this class.
  The plot method displays the entire array of functions.
  The method \code{\link{[.fasp}} selects a sub-array using the natural
  indices \code{i,j}. 

  An object of class \code{"fasp"} is a list containing at least the
  following components:
  
  \item{fns}{
    A list of data frames, each representing one of the functions.
  }
  \item{which}{
    A matrix representing the spatial arrangement of the
    functions. If \code{which[i,j] = k}
    then the function represented by \code{fns[[k]]} should be plotted
    in the panel at position \eqn{(i,j)}. If \code{which[i,j] = NA}
    then nothing is plotted in that position.
  }
  \item{titles}{
    A list of character strings, providing suitable plotting titles
    for the functions.
    }
  \item{default.formulae}{
    A list of default formulae for plotting each of the functions.
  }
  \item{title}{
    A character string, giving a default title for the array
    when it is plotted.
  }

}
\seealso{
  \code{\link{alltypes}},
  \code{\link{allstats}},
  \code{\link{plot.fasp}},
  \code{\link{[.fasp}}
}
\examples{
  require(spatstat)
  # unmarked point pattern
  data(swedishpines)
  a <- allstats(swedishpines,dataname="Swedish Pines")
  a
  plot(a)
  plot(a[1,])

  # multitype point pattern
  data(amacrine)
  a <- alltypes(amacrine, "G")
  plot(a)

  # select the row corresponding to cells of type "on"
  b <- a["on", ]
  plot(b)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 
