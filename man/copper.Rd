\name{copper}
\alias{copper}
\title{
  Berman-Huntington points and lines data
}
\description{
  These data come from an intensive geological survey of
  a 70 x 158 km region in central Queensland, Australia.
  They consist of 67 points representing copper ore deposits,
  and 146 line segments representing geological `lineaments'.
  Lineaments are linear features, visible on a satellite image,
  that are believed to consist largely of geological faults (Berman, 1986,
  p. 55).
  It would be of great interest to predict the occurrence of copper deposits
  from the lineament pattern, since the latter can easily be observed on
  satellite images. 

  These data were introduced and analysed by Berman (1986).
  They have also been studied by Berman and Turner (1992),
  Baddeley and Turner (2000) and Foxall and Baddeley (2002).
  
  Many analyses have been performed on the southern half of the data only.
  This subset is also provided.
} 
\format{
  \code{copper} is a list with the following entries:
  \describe{
    \item{points}{an object of class \code{"ppp"}
      representing the full point pattern of copper deposits.
      See \code{\link{ppp.object}} for details of the format.
    }
    \item{lines}{the coordinates of the lineaments in the full dataset,
      given as a data frame with 4 columns (x1, y1, x2, y2).
    }
    \item{SouthWindow}{the window delineating the southern half of
      the study region. An object of class \code{"owin"}.
    }
    \item{SouthPoints}{the point pattern of copper deposits in the
      southern half of the study region. An object of class
      \code{"ppp"}.
    }
    \item{SouthLines}{the coordinates of the lineaments in the
      southern half of the study region. A data frame with 4 columns
      in format (x1, y1, x2, y2).
    }
    \item{SouthDistance}{A function with no arguments. The return value
      is a pixel image with greyscale value equal to the
      distance to the nearest lineament. Computed only inside the
      southern half of the dataset. An object of class \code{"im"}.
    }
  }
}
\usage{data(copper)}
\examples{

  data(copper)

  # Plot full dataset

  plot(copper$points)
  cl <- copper$lines
  segments(cl[,1], cl[,2], cl[,3], cl[,4])

  # Plot southern half of data
  plot(copper$SouthPoints)
  cl <- copper$SouthLines
  segments(cl[,1], cl[,2], cl[,3], cl[,4])

  \dontrun{
    Z <- copper$SouthDistance()
    plot(Z)
    X <- copper$SouthPoints
    ppm(X, ~D, covariates=list(D=Z))
  }
}
\source{
  Dr J. Huntington.
  Coordinates kindly provided by Dr. Mark Berman, CSIRO, Sydney,
  Australia.
}
\references{
Baddeley, A. and Turner, R. (2000)
Practical maximum pseudolikelihood for spatial point patterns.
\emph{Australian and New Zealand Journal of Statistics}
\bold{42}, 283--322.
 
Berman, M. (1986).
Testing for spatial association between a point process and another
  stochastic process.
\emph{Applied Statistics} \bold{35}, 54--62.

Berman, M. and Turner, T.R. (1992)
Approximating point process likelihoods with GLIM.
\emph{Applied Statistics} \bold{41}, 31--38.
 
Foxall, R. and Baddeley, A. (2002)
Nonparametric measures of association between a
spatial point process and a random set, with
geological applications. \emph{Applied Statistics} \bold{51}, 165--182.

}
\keyword{datasets}
