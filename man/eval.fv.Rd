\name{eval.fv}
\alias{eval.fv}
\title{Evaluate Expression Involving Functions}
\description{
  Evaluates any expression involving one or more function value (fv) objects,
  and returns another object of the same kind.
}
\usage{
  eval.fv(expr, envir)
}
\arguments{
  \item{expr}{An expression.}
  \item{envir}{Optional. The environment in which to evaluate the expression.}
}
\details{
  This is a wrapper to make it easier to perform
  pointwise calculations with the summary functions
  used in spatial statistics.

  An object of class \code{"fv"} is essentially a data frame
  containing several different statistical estimates of the same
  function. Such objects are returned by \code{\link{Kest}} and its
  relatives.

  For example, suppose \code{X} is an object of class \code{"fv"}
  containing several different estimates of the Ripley's K function \eqn{K(r)},
  evaluated at a sequence of values of \eqn{r}.
  Then \code{eval.fv(X+3)} effectively adds 3 to 
  each function estimate in \code{X}, and returns
  the resulting object. 

  Suppose \code{X} and \code{Y} are two objects of class \code{"fv"}
  which are compatible (in particular they have the same vector
  of \eqn{r} values). Then 
  \code{eval.im(X + Y)} will add the corresponding function values in
  \code{X} and \code{Y}, and return the resulting function.

  In general, \code{expr} can be any expression involving
  (a) the \emph{names} of objects of class \code{"fv"}, (b) scalar
  constants, and (c) functions which are vectorised.
  See the Examples.

  First \code{eval.fv} determines which of the \emph{variable names}
  in the expression \code{expr} refer to objects of class \code{"fv"}.
  Each such name is replaced by a vector containing the function values.
  The expression is then evaluated. The result should be a vector;
  it is taken as the new vector of function values.

  The expression \code{expr} must be vectorised.
  There must be at least one object of class \code{"fv"} in the expression.
  All such objects must be compatible.
}
\value{
  Another object of class \code{"fv"}.
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{Kest}}
}
\examples{
  # manipulating the K function
  X <- rpoispp(42)
  Ks <- Kest(X)

  eval.fv(Ks + 3)
  Ls <- eval.fv(sqrt(Ks/pi))

  # manipulating two K functions
  Y <- rpoispp(20)
  Kr <- Kest(Y)
  
  Kdif <- eval.fv(Ks - Kr)
  Z <- eval.fv(sqrt(Ks/pi) - sqrt(Kr/pi))
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
