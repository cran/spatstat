\name{model.depends}
\alias{model.depends}
\alias{model.is.additive}
\title{
  Identify Covariates Involved in each Model Term
}
\description{
  Given a fitted model (of any kind), identify which of the covariates
  is involved in each term of the model.
}
\usage{
model.depends(object)
model.is.additive(object)
}
\arguments{
  \item{object}{
    A fitted model of any kind.
  }
}
\details{
  The \code{object} can be a fitted model of any kind,
  including models of the classes \code{\link{lm}}, \code{\link{glm}}
  and \code{\link{ppm}}.

  To be precise,
  \code{object} must belong to a class for which there are methods
  for \code{\link{formula}}, \code{\link{terms}}
  and \code{\link{model.matrix}}.
  
  The command \code{model.depends} determines the relationship between
  the original covariates (the data supplied when \code{object} was
  fitted) and the canonical covariates (the columns of the design matrix).
  It returns a logical matrix, with one row for each canonical
  covariate, and one column for each of the original covariates,
  with the \code{i,j} entry equal to \code{TRUE} if the
  \code{i}th canonical covariate depends on the \code{j}th
  original covariate.

  The command \code{model.is.additive} determines whether the model
  is additive, in the sense that there is no canonical covariate that
  depends on two or more original covariates. It returns a logical value.
}
\value{
  Logical value or matrix.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{model.matrix}}
}
\examples{
   x <- 1:10
   y <- 3*x + 2
   z <- rep(c(-1,1), 5)
   fit <- lm(y ~ poly(x,2) + sin(z))
   model.depends(fit)
   model.is.additive(fit)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}

