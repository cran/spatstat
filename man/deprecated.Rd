\name{spatstat-deprecated}
\alias{conspire}
\alias{eval.hyper}
\alias{ksmooth.ppp}
\alias{ks.test.ppm}
\alias{mpl}
\alias{rtoro}
\title{Deprecated spatstat functions}
\description{
  Deprecated spatstat functions.
}
\usage{
conspire(\dots)
ksmooth.ppp(x, sigma, \dots, edge=TRUE)
ks.test.ppm(\dots)
mpl(Q, trend, interaction, data, correction, rbord, use.gam)
rtoro(X, which=NULL, radius=NULL, width=NULL, height=NULL) 
eval.hyper(e, h, simplify=TRUE, ee=NULL)
}

\details{
  \code{conspire} has been replaced by \code{\link{plot.fv}}.

  \code{ksmooth.ppp} has been replaced by \code{\link{density.ppp}}.

  \code{mpl} has been replaced by \code{\link{ppm}}.

  \code{ks.test.ppm} has been replaced by \code{\link{kstest.ppm}}.
  
  \code{rtoro} has been replaced by \code{\link{rshift}}.

  \code{eval.hyper} has been replaced by \code{\link{with.hyperframe}}.
}
\keyword{internal}
