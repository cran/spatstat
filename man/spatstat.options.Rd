\name{spatstat.options}
\alias{spatstat.options}
\title{Internal Options in Spatstat Package}
\description{
  Allows the user to examine and reset the values
  of global parameters which control actions in the
  \code{spatstat} package.
}
\usage{
  spatstat.options(...)
}
\arguments{
  \item{\dots}{
    Either empty,
    or a succession of parameter names in quotes,
    or a succession of \code{name=value} pairs.
    See below for the parameter names.
  }
}
\value{
  A list of parameters and their values.
}
\details{
  This function is analogous to \code{\link{options}}
  and calls the utility function \code{\link{check.options}}.

  If no arguments are given, 
  the current values of all parameters are returned.

  If a succession of parameter names is given,
  the current values of these parameters are returned.
  
  If \code{name=value} pairs are given, the named parameters
  are reset to the given values, and the \bold{previous} values of
  these parameters are returned.

  The parameters are:
  \describe{
    \item{npixel}{
      Controls the number of pixels in the discrete approximation
      of an irregular window by a binary pixel image as created by
      \code{\link{as.mask}}. Either an integer or a pair of integers
      giving the number of pixels in the \code{x} and \code{y}
      directions.
    }
    \item{maxedgewt}{
      Edge correction weights will be trimmed
      so as not to exceed this value.
      This applies to the weights computed by
      \code{\link{edge.Trans}} or \code{\link{edge.Ripley}}
      and used in \code{\link{Kest}} and its relatives.
    }
    \item{par.binary}{
      List of arguments to be passed to the function \code{\link{image}}
      when displaying a binary image mask (in \code{\link{plot.owin}}
      or \code{\link{plot.ppp}}).
      Typically used to reset the colours of foreground and background.
    }
    \item{par.persp}{
      List of arguments to be passed to the function \code{\link{persp}}
      when displaying a real-valued image, such as the fitted surfaces
      in \code{\link{plot.ppm}}.
    }
  }

}
\seealso{
  \code{\link{options}}
}
\examples{
  \dontrun{library(spatstat)}

  spatstat.options(npixel=150)
  spatstat.options(npixel=c(100,200))

  spatstat.options(par.binary=list(col=grey(c(0.5,1))))

  spatstat.options(par.persp=list(theta=-30,phi=40,d=4))
  # see help(persp.default) for other options

}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
