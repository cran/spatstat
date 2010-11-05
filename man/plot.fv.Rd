\name{plot.fv}
\alias{plot.fv}
\title{Plot Function Values}
\description{
  Plot method for the class \code{"fv"}.
}
\usage{
 \method{plot}{fv}(x, fmla, \dots, subset=NULL, lty=NULL, col=NULL, lwd=NULL,
           xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, ylim.covers=NULL,
           legend=TRUE, legendpos="topleft", legendmath=FALSE,
           shade=NULL, shadecol="grey")
}
\arguments{
  \item{x}{
    An object of class \code{"fv"}, containing the variables to be plotted
    or variables from which the plotting coordinates can be computed.
  }
  \item{fmla}{
    an R language formula 
    determining which variables or expressions are plotted.
    Either a formula object, or a string that can be parsed as a formula.
  }
  \item{subset}{
    (optional) subset of rows of the data frame that will be plotted.
  }
  \item{lty}{
    (optional) numeric vector of values of the graphical parameter
    \code{lty} controlling the line style of each plot.
  }
  \item{col}{
    (optional) numeric vector of values of the graphical parameter
    \code{col} controlling the colour of each plot.
  }
  \item{lwd}{
    (optional) numeric vector of values of the graphical parameter
    \code{lwd} controlling the line width of each plot.
  }
  \item{xlim}{
    (optional) range of x axis
  }
  \item{ylim}{
    (optional) range of y axis
  }
  \item{xlab}{
    (optional) label for x axis
  }
  \item{ylab}{
    (optional) label for y axis
  }
  \item{\dots}{
    Extra arguments passed to \code{plot.default}.
  }
  \item{ylim.covers}{
    Optional vector of \eqn{y} values that must be included in the
    \eqn{y} axis. For example \code{ylim.covers=0} will ensure that the
    \eqn{y} axis includes the origin.
  }
  \item{legend}{
    Logical flag or \code{NULL}. If \code{legend=TRUE}, the algorithm
    plots a legend in the top left corner of the plot,
    explaining the meaning of the different line types and colours.
  }
  \item{legendpos}{
    The position of the legend. Either a character string keyword
    (see \code{\link[graphics]{legend}} for keyword options)
    or a pair of coordinates in the format \code{list(x,y)}.
  }
  \item{legendmath}{
    Logical. If \code{TRUE}, the legend will display the
    mathematical notation for each curve. If \code{FALSE}, the legend text
    is the identifier (column name) for each curve.
  }
  \item{shade}{
    An index that identifies two columns of \code{x}. When the corresponding
    curves are plotted, the region between the curves will be 
    shaded in light grey. 
    Often used for displaying simulation envelopes, by setting
    \code{shade=c("hi", "lo")}.
  }
  \item{shadecol}{
    The colour to be used in the \code{shade} plot.
    A character string or an integer specifying a colour.
  }
}
\value{
  Either \code{NULL}, or a data frame giving the meaning of the
  different line types and colours.
}
\details{
  This is the \code{plot} method for the class \code{"fv"}.

  The use of the argument \code{fmla} is like \code{plot.formula}, but offers
  some extra functionality.
  
  The left and right hand sides of \code{fmla}
  are evaluated in the data frame \code{x},
  and the results are plotted against each other
  (the left side on the \eqn{y} axis 
  against the right side on the \eqn{x} axis). 
  Both left and right sides may be variables in the data frame
  or expressions in these variables.

  Multiple curves may be specified by a single formula
  of the form 
  \code{cbind(y1,y2,\dots,yn) ~ x}, where \code{x,y1,y2,\dots,yn} are
  expressions involving the variables in the data frame.
  Each of the variables \code{y1,y2,\dots,yn} in turn will be plotted
  against \code{x}. 
  See the examples.

  A convenient abbreviation is the symbol \code{.} which
  can be used in the formula to represent all
  the variables in the data frame (other than the function argument
  itself). See the examples.

  The value returned by this plot function indicates the
  meaning of the line types and colours in the plot. It can be used
  to make a suitable legend for the plot if you want to do this
  by hand. See the examples.

  The argument \code{shade} can be used to display critical bands
  or confidence intervals. If it is not \code{NULL}, then it should be
  a subset index for the columns of \code{x}, that identifies exactly
  2 columns. When the corresponding curves are plotted, the region
  between the curves will be shaded in light grey. See the Examples.
}
\examples{
   data(cells)
   K <- Kest(cells)
   # K is an object of class "fv"

   plot(K, iso ~ r)                # plots iso against r

   plot(K, sqrt(iso/pi) ~ r)   # plots sqrt(iso/r)  against r

   plot(K, cbind(iso,theo) ~ r)   # plots iso against r  AND theo against r

   plot(K, .  ~ r)            # plots all available estimates of K against r

   plot(K, sqrt(./pi) ~ r)   # plots all estimates of L-function
                             # L(r) = sqrt(K(r)/pi)

   plot(K, cbind(iso,theo) ~ r, col=c(2,3))
                                   # plots iso against r  in colour 2
                                   # and theo against r in colour 3

   plot(K, iso ~ r, subset=quote(r < 0.2))
                                   # plots iso against r for r < 10

   # making a legend by hand
   v <- plot(K, . ~ r, ylab="K(r)", legend=FALSE)
   legend(0.05, 0.15, legend=v$meaning, lty=v$lty, col=v$col)

   # significance bands
   KE <- envelope(cells, Kest, nsim=19)
   plot(KE, shade=c("hi", "lo"))
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{Kest}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
 
