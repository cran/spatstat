\name{ganglia}
\alias{ganglia}
\title{Cat Retinal Ganglia Data} 
\description{
  Point pattern of retinal ganglion cells identified as `on' or `off'.
  A marked point pattern.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of cell \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of cell \cr
    \code{marks} \tab factor with levels \code{off} and \code{on} \cr
                 \tab indicating ``off'' and ``on'' cells
  }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(ganglia)}
\source{W\"assle et al (1981), data supplied by Peter Diggle}
\section{Notes}{
The data represent a pattern of beta-type ganglion cells
in the retina of a cat recorded by W\"assle et al. (1981).
The window is a rectangle 
\eqn{(0,1) \times (0, .7533)}{(0,1) x (0, .7533)} mm approximately.

Beta cells are associated
with the resolution of fine detail in the cat's visual system.
They can be classified anatomically as ``on'' or ``off''.
Statistical independence of the ``on''- and ``off''-components would strengthen
the claim that there are two separate channels for 
``brightness'' and ``darkness'' as postulated by Hering in 1874.
For details see Waessle et al (1981).

These data have been analysed by Diggle (1986) and
Van Lieshout and Baddeley (1999).
}
\references{
Diggle, P. J. (1986).
Displaced amacrine cells in the retina of a
rabbit : analysis of a bivariate spatial point pattern. 
\emph{J. Neurosci. Meth.} \bold{18}, 115--125.
 
Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
Indices of dependence between types in multivariate point patterns.
\emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

W\"assle, H., Boycott, B. B. & Illing, R.-B. (1981).
Morphology and mosaic of on- and off-beta cells in the cat retina and
some functional considerations.
\emph{Proc. Roy. Soc. London Ser. B} \bold{212}, 177--195.
}
\keyword{datasets}
