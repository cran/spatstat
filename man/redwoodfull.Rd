\name{redwoodfull}
\alias{redwoodfull}
\alias{redwoodfull.diag}
\alias{redwoodI}
\alias{redwoodII}
\title{
  California Redwoods Point Pattern (Entire Dataset)
}
\description{
  Locations of 195 seedlings and saplings 
  of California redwood trees.
 
  The data \code{redwoodfull}
  represent the locations of 195 seedlings and saplings 
  of California redwood trees in a square sampling region.  
  They were described and analysed by Strauss (1975).

  Strauss divided the sampling region into two subregions I and II
  demarcated by a diagonal line across the region. The dataset
  \code{redwoodfull.diag} specifies this diagonal line,
  and the datasets \code{redwoodI} and \code{redwoodII}
  represent the point patterns in these two subregions. 

  Yet another subset of the full dataset, consisting of 62 points
  in a square subregion, was extracted by Ripley (1977). It
  has been re-analysed many times, and is the dataset usually known as
  "the redwood data" in the spatial statistics literature.
  That dataset is supplied in the \code{spatstat}
  library as \code{\link{redwood}}.

} 
\format{
  The dataset \code{redwoodfull} is an object of class \code{"ppp"}
  representing the point pattern of tree locations.
  The window has been rescaled to the unit square.
  See \code{\link{ppp.object}} for details of the format of a
  point pattern object.

  The datasets \code{redwoodI} and \code{redwoodII} are also objects of
  class \code{"ppp"} representing the point patterns in regions I and II
  respectively, on the same scale as \code{redwoodfull}.

  The dataset \code{redwoodfull.diag} is a list with components
  \code{x} and \code{y} giving the \eqn{x,y} coordinates of the
  two endpoints of the diagonal line separating regions I and II.
  Its main use is in plotting the data.
}
\usage{data(redwoodfull)
       plot(redwoodfull)
       lines(redwood.diag)
       plot(redwoodI)
       plot(redwoodII)
}
\source{The plot of the data published by Strauss (1975)
  was scanned and digitised by Sandra Pereira, University of
  Western Australia, 2002.
}
\seealso{
  \code{\link{redwood}}
}
\references{
  Diggle, P.J. (1983)
  \emph{Statistical analysis of spatial point patterns}.
  Academic Press.

  Ripley, B.D. (1977)
  Modelling spatial patterns (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{39}, 172--212.

  Strauss, D.J. (1975)
  A model for clustering.
  \emph{Biometrika} \bold{63}, 467--475.
}
\keyword{datasets}
