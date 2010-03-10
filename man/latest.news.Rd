\name{latest.news}
\alias{latest.news}
\title{
  Print News About Latest Version of Package
}
\description{
  Prints the news documentation for the current version of \code{spatstat}
  or another specified package.
}
\usage{
latest.news(package = "spatstat")
}
\arguments{
  \item{package}{
    Name of package for which the news file
  }
}
\details{
  By default, this function prints the news documentation about changes in the
  current installed version of the \pkg{spatstat} package.

  If \code{package} is given, then the function reads the
  news for the specified package from its \code{NEWS} file (if it has one) 
  and prints only the entries that refer to the current version
  of the package.

  To see the news for all previous versions as well as the current
  version, use the \R utility \code{\link[utils]{news}}. See the Examples.
}
\value{
  Null.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link[utils]{news}}
}
\examples{
  # current news
  latest.news()
  # all news
  news(package="spatstat")
}
\keyword{documentation}