\name{eval.hyper}
\alias{eval.hyper}
\title{Evaluate an Expression in Each Row of a Hyperframe}
\description{
  An expression, involving the names of columns in a hyperframe,
  is evaluated separately for each row of the hyperframe.
}
\usage{
eval.hyper(e, h, simplify = TRUE, ee = NULL)
}
\arguments{
  \item{e}{An \R language expression to be evaluated.}
  \item{h}{A hyperframe (object of class \code{"hyperframe"})
    containing data.
  }
  \item{simplify}{
    Logical. If \code{TRUE}, the return value
    will be simplified to a vector whenever possible.
  }
  \item{ee}{
    Alternative form of \code{e}. See Details.
  }
}
\details{
  This function evaluates the expression \code{e} in each row
  of the hyperframe \code{h}.
  
  The argument \code{e} should be an \R language expression
  in which each variable name is either the name of a column in the
  hyperframe \code{h}, or the name of an object in the global environment.
  
  For each row of \code{h}, the expression will be evaluated
  so that variables which are column names of \code{h} are
  interpreted as the entries for those columns in the current row.

  For example, if \code{h} has columns called \code{A} and \code{B},
  then \code{eval.hyper(A != B, h)} inspects each row of \code{h} in turn,
  tests whether the entries in columns \code{A} and \code{B} are
  equal, and returns the \eqn{n} logical values.


}
\value{
  Normally a list of length
  \eqn{n} (where \eqn{n} is the number of rows) containing the results
  of evaluating the expression for each row. 
  If \code{simplify=TRUE} and each result is a single atomic value,
  then the result is a vector or factor 
  containing the same values.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\seealso{
  \code{\link{hyperframe}},
  \code{\link{plot.hyperframe}}
}
\examples{
  # generate Poisson point patterns with intensities 10 to 100
   H <- hyperframe(L=seq(10,100, by=10))
   X <- eval.hyper(rpoispp(L),H)
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
