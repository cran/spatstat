\name{summary.ppm}
\alias{summary.ppm}
\alias{print.summary.ppm}
\title{Summarizing a Fitted Point Process Model}
\description{
  \code{summary} method for class \code{"ppm"}.
}
\usage{
  summary.ppm(object, \dots, quick=FALSE)
  print.summary.ppm(x, \dots)
}
\arguments{
  \item{object}{A fitted point process model.}
  \item{\dots}{Ignored.}
  \item{quick}{Logical flag controlling the scope of the summary.}
  \item{x}{Object of class \code{"summary.ppm"} as returned by
    \code{summary.ppm}.
  }
}
\details{
  This is a method for the generic \code{\link{summary}}
  for the class \code{"ppm"}. An object of class \code{"ppm"}
  describes a fitted point process model. See \code{\link{ppm.object}})
  for details of this class. 
  
  \code{summary.ppm} extracts information about the
  type of model that has been fitted, the data to which the model was
  fitted, and the values of the fitted coefficients.
  (If \code{quick=TRUE} then only the information about the type
  of model is extracted.)

  \code{print.summary.ppm} prints this information in a
  comprehensible format.

  In normal usage, \code{print.summary.ppm} is invoked implicitly
  when the user calls \code{summary.ppm} without assigning its value
  to anything. See the examples.
}
\value{
  \code{summary.ppm} returns an object of class \code{"summary.ppm"},
  while \code{print.summary.ppm} returns \code{NULL}.
}
\examples{
  # invent some data
  X <- rpoispp(42)
  # fit a model to it
  fit <- ppm(X, ~x, Strauss(r=0.1))
  # summarize the fitted model
  summary(fit)
  # `quick' option
  summary(fit, quick=TRUE)

  # save the full summary
  s <- summary(fit)
  # print it
  print(s)
  s
  # extract stuff
  names(s)
  s$args$correction
  s$name
  s$trend$value

  # multitype pattern
  data(demopat)
  fit <- ppm(demopat, ~marks, Poisson())
  summary(fit)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
