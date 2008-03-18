\name{progressreport}
\alias{progressreport}
\title{Print Progress Reports}
\description{
  Prints Progress Reports during a loop or iterative calculation.
}
\usage{
progressreport(i, n, every = max(1, ceiling(n/100)), nperline = min(charsperline, every * ceiling(charsperline/(every + 3))), charsperline = 60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
    Integer. The current iteration number (from 1 to \code{n}).
  }
  \item{n}{
    Integer. The (maximum) number of iterations to be computed.
  }
  \item{every}{
    Optional integer. The number of iterations between successive
    reports.
  }
  \item{nperline}{
    Optional integer. The maximum number of reports to be printed per
    line of output.
  }
  \item{charsperline}{
    Optional integer. The number of characters in a line of output.
  }
}
\details{
  This is a convenient function for reporting progress
  during an iterative sequence of calculations
  or a suite of simulations. 
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
\examples{
  for(i in 1:40) progressreport(i, 40)
}
\keyword{print}
