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
  split.ppp(x, f = x$marks, drop=FALSE, un=NULL, \dots)
}
\usage{
  split.ppp(x, f = x$marks, drop=FALSE, un=missing(f), \dots)
  split(x, f, drop, un) <- value
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
  }
  \item{f}{
    Factor determining the grouping. 
  }
  \item{drop}{
	Logical. Determines whether empty groups will
	be deleted.
  }
  \item{un}{
    Logical. Determines whether subpatterns will be unmarked
    (i.e. whether marks will be removed	from the points in each
    subpattern). 
  }
  \item{\dots}{
    Other arguments are ignored.
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

  Some of the sub-patterns created by the split
  may be empty. If \code{drop=TRUE}, then empty sub-patterns will
  be deleted from the list. If \code{drop=FALSE} then they are retained.

  The argument \code{un} determines how to handle marks 
  in the case where \code{x} is a marked point pattern.
  If \code{un=TRUE} then the marks of the 
  points will be discarded when they are split into groups,
  while if \code{un=FALSE} then the marks will be retained.
  
  The result of \code{split.ppp} has class \code{"splitppp"}
  and can be plotted using \code{\link{plot.splitppp}}.

  The assignment function \code{split<-.ppp} 
  updates the point pattern \code{x} so that
  it satisfies \code{split(x, f, drop, un) = value}. The argument \code{value}
  is expected to be a list of point patterns, one for each level of
  \code{f}. These point patterns are expected to be compatible in the sense
  that they all have the same window, and either they are all unmarked
  or they all have marks of the same kind. 

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

# plot them
 plot(split(amacrine))

# the following are equivalent:
 amon <- split(amacrine)$on
 amon <- unmark(amacrine[amacrine$marks == "on"])
   
# the following are equivalent:
 amon <- split(amacrine, un=FALSE)$on
 amon <- amacrine[amacrine$marks == "on"]
   
# Scramble the locations of the 'on' cells
 u <- split(amacrine)
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
