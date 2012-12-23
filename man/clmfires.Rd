\name{clmfires}
\alias{clmfires}
\alias{clmcov100}
\alias{clmcov200}
\docType{data}
\title{
  Castilla-La Mancha Forest Fires
}
\description{
  This dataset is a record of forest fires
  in the Castilla-La Mancha region of Spain
  between 1997 and 2007. 
  
  The dataset \code{clmfires} is a point pattern
  (object of class \code{"ppp"}) containing the spatial coordinates
  of each fire, with marks containing information about each fire.
  There are 11 columns of marks:
  \tabular{ll}{
    \code{cause} \tab cause of fire (see below) \cr
    \code{burnt.area} \tab total area burned, in hectares \cr
    \code{year} \tab calendar year of fire (factor with levels 1997-2007) \cr
    \code{month} \tab calendar month of fire (factor with levels 1-12) \cr
    \code{day} \tab calendar day of fire (integer 1-31) \cr
    \code{time} \tab number of days elapsed since 1 January 1997 \cr
    \code{elevation} \tab terrain elevation in metres above sea level\cr
    \code{orientation} \tab direction of upward slope in degrees (0-360) \cr
    \code{slope} \tab slope of terrain \cr
    \code{slope.discretized} \tab discretized value of \code{slope} \cr
    \code{land.use} \tab type of land cover (see below)
  }
  The \code{cause} of the fire is a factor with the levels
  \code{lightning}, \code{accident} (for accidents or negligence),
  \code{intentional} (for intentionally started fires) and 
  \code{unknown} (for other causes or unknown cause).

  The \code{land.use} variable is a factor with levels
  \code{urban}, \code{farm} (for farms or orchards),
  \code{meadow}, \code{denseforest} (for dense forest),
  \code{conifer} (for conifer forest or plantation),
  \code{mixedforest}, \code{grassland}, \code{bush}, \code{scrub}
  and \code{artifgreen} for artificial greens such as golf courses.
  
  Accompanying this point pattern, there are two datasets
  \code{clmcov100} and \code{clmcov200} containing covariate information
  for the entire Castilla-La Mancha region. Each of these two datasets
  is a list of four images (objects of class \code{"im"})
  named \code{elevation}, \code{orientation}, \code{slope} and
  \code{landuse}. These images provide values for the four covariates
  at every location in the study area. The images in  \code{clmcov100}
  are 100 by 100 pixels in size, while those in  \code{clmcov200} are
  200 by 200 pixels.
  For easy handling, \code{clmcov100} and \code{clmcov200}
  also belong to the class \code{"listof"} so that they can be
  plotted and printed immediately.
} 
\format{
  \code{clmfires} is a marked point pattern (object of class \code{"ppp"}). 
  See \code{\link[spatstat]{ppp.object}}.

  \code{clmcov100} and \code{clmcov200} are lists of pixel images
  (objects of class \code{"im"}).
}
\usage{data(clmfires)}
\examples{
plot(clmfires, which.marks="cause", cols=2:5, cex=0.25)
plot(clmcov100)
}
\source{
  Jorge Mateu.
}
\keyword{datasets}
\keyword{spatial}
