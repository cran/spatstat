\name{split.ppp}
\alias{split.ppp}
\alias{split<-.ppp}
\alias{split}
\alias{split<-}
\title{Divide Point Pattern into Sub-patterns}
\description{
  Divides a point pattern into several sub-patterns,
  according to their marks, or according to any user-specified grouping.
}
\synopsis{
  split.ppp(x, f = x$marks)
}
\usage{
  split.ppp(x, f = x$marks)
  split(x, f) <- value
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
  }
  \item{f}{
    Factor determining the grouping. 
  }
  \item{value}{
    List of point patterns.
  }
} 
\value{
  The value of \code{split.ppp} is a list of point patterns.
  The components of the list are named by the levels of \code{f}.

  The assignment form \code{split<-.ppp} returns the updated
  point pattern \code{x}. 
}
\details{
  The function \code{split.ppp}
  divides up the points of the point pattern \code{x}
  into several sub-patterns according to the levels of the factor
  \code{f}. The result is a list of point patterns, one for each
  level of \code{f}.

  If \code{f} is present, it must be a factor, and its length
  must equal the number of points in \code{x}. The levels of \code{f}
  determine the destination of each point in \code{x}.
  The \code{i}th point of \code{x} will be placed in the sub-pattern
  \code{split.ppp(x)$l} where \code{l = f[i]}.

  If \code{f} is missing, then \code{x} must be a multitype point pattern
  (a marked point pattern whose marks vector is a factor).
  Then the effect is that the points of each type
  are separated into different point patterns.

  The result of \code{split.ppp} has class \code{"splitppp"}
  and can be plotted using \code{\link{plot.splitppp}}.

  The assignment function \code{split<-.ppp} 
  updates the point pattern \code{x} so that
  it satisfies \code{split(x, f) = value}. The argument \code{value}
  is expected to be a list of point patterns, one for each level of
  \code{f}. 

  Splitting can also be undone by the function \code{\link{superimpose}}.
}
\seealso{
  \code{\link{cut.ppp}},
  \code{\link{plot.splitppp}},
  \code{\link{superimpose}},
  \code{\link{ppp.object}}
}
\examples{

# Multitype point pattern: separate into types
 data(amacrine)
 u <- split(amacrine)

# the following are equivalent:
 amon <- amacrine[amacrine$marks == "on"]
 amon <- split(amacrine)$on
   
# plot them
 plot(split(amacrine))

# Scramble the 'on' cells
 u$on <- runifpoint(amon$n, amon$window)
 split(amacrine) <- u

# Point pattern with continuous marks
 data(longleaf)
 \testonly{
	# smaller dataset
	longleaf <- longleaf[seq(1, longleaf$n, by=80)]
 }
 # cut the range of tree diameters into three intervals
 long3 <- cut.ppp(longleaf, 3)
 # now split them
 long3split <- split(long3)

# Unmarked point pattern
  data(swedishpines)
# cut & split according to nearest neighbour distance
  f <- cut(nndist(swedishpines), 3)
  u <- split(swedishpines, f)
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
