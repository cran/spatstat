\name{plot.ppp}
\alias{plot.ppp}
\title{plot a Spatial Point Pattern}
\description{
  Plot a two-dimensional spatial point pattern
}
\usage{
 plot.ppp(x, main, \dots, chars, cols, use.marks=TRUE, add=FALSE, maxsize)
}
\arguments{
  \item{x}{
    The spatial point pattern to be plotted.
    An object of class \code{"ppp"},
    or data which can be converted into 
    this format by \code{\link{as.ppp}()}.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting functions
    \code{\link{plot.default}} (see \code{\link{par}} for possible
    options) or \code{\link{symbols}}.
  }
  \item{chars}{
    the plotting characters which should be used to represent
    the points of different types, in the case of a multitype
    point pattern. 
  }
  \item{cols}{
    the colours which should be used to plot the points of
    different types, in the case of a multitype point pattern.
  }
  \item{use.marks}{
    logical flag; if \code{TRUE}, plot points using a different
    plotting symbol for each mark;
    if \code{FALSE}, only the locations of the points will be plotted,
    using \code{\link{points}()}.
  }
  \item{add}{
    logical flag; if \code{TRUE},
    just the points are plotted, over the existing plot.
    A new plot is not created, and
    the window is not plotted.
  }
  \item{maxsize}{
    the maximum radius of the circles plotted.
    If \code{x} is a marked point pattern with 
    numerical marks, the marks are plotted as circles with
    maximum radius \code{maxsize} 
    (in the same units as the coordinates of
    the point pattern). There is a sensible default.
  }
}
\value{
  \code{NULL}, or a vector giving the correspondence between 
  mark values and plotting characters.
}
\details{
  This is the \code{plot} method for 
  point pattern datasets (of class \code{"ppp"}, see \code{\link{ppp.object}}).

  First the observation window \code{x$window} is plotted 
  using \code{\link{plot.owin}}.
  Then the locations of the points are plotted. 
 
  If the points do not have marks (i.e. \code{x$marks} is absent)
  or if \code{use.marks = FALSE},
  then the locations of the points will be plotted 
  using a single plot character. This character can be changed
  by the extra argument \code{pch};  see the examples.

  For a marked point pattern (where \code{x$marks} is not null)
  each point will be represented by a graphical symbol plotted at the
  location of the point.
  The graphical representation of the marks depends
  on their storage mode:
  \describe{
    \item{multitype point pattern:}{
      If \code{x$marks} is a factor, then 
      each level of the factor is 
      represented by a different plot character.
    }
    \item{continuous marks:}{
      If \code{x$marks} is a numeric vector,
      the marks are rescaled to the unit interval and
      each point is represented by a circle 
      with radius proportional to the rescaled mark
      (if the value is positive) or a square with side proportional
      to the absolute value of the rescaled mark (if the value is negative).
    }
    \item{otherwise:}{
      Each possible mark will be represented by a
      different plotting character. The default is to
      represent the \eqn{i}th smallest mark value by
      \code{points(..., pch=i)}.
    }
  }
  If the arguments \code{chars} and/or \code{cols} are supplied,
  they should be vectors of the same length
  as \code{sort(unique(x$marks))}, and then the \eqn{i}th smallest mark
  will be plotted using character \code{chars[i]} and colour
  \code{cols[i]}.

  The argument \code{cols} is incompatible with
  the generic plot argument \code{col}. If both are given,
  an error is generated and no plot is produced.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{plot}},
  \code{\link{par}},
  \code{\link{points}},
  \code{\link{plot.owin}},
  \code{\link{symbols}}
}
\examples{
   \dontrun{
   data(cells)
   plot(cells)

   # multitype 
   data(lansing)
   plot(lansing)

   # marked by a real number
   data(longleaf)
   plot(longleaf)

   # just plot the points
   plot(longleaf, use.marks=FALSE)
   plot(unmark(longleaf)) # equivalent

   # controlling COLOURS
   plot(cells, col="blue")
   plot(lansing, cols=c("black", "yellow", "green", "pink",
                        "blue","red","white")
   plot(longleaf, fg="blue")

   # use selected plot characters
   lma <- sort(unique(longleaf$marks))
   adult <- ifelse(lma >= 30, "o", "+")
   plot(longleaf, chars = adult)

   # make the plotting symbols larger (for publication at reduced scale)
   plot(cells, cex=2)

  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
