\name{allstats}
\alias{allstats}
\title{Calculate four standard summary functions of a point pattern.}
\description{
  Calculates the \eqn{F}, \eqn{G}, \eqn{J}, and \eqn{K}
  summary functions for an unmarked point pattern.
  Returns them as a function array (of class \code{"fasp"}, see
  \code{\link{fasp.object}}).
}
\usage{
  allstats(pp, dataname=NULL, verb=F)
}
\arguments{
  \item{pp}{The observed point pattern, for which summary function
    estimates are required.  An object of class \code{"ppp"}.
    It must not be marked.
  }
  \item{dataname}{A character string giving an optional (alternative)
    name for the point pattern.
  }
  \item{verb}{A logical value meaning ``verbose''. If \code{TRUE},
    progress reports are printed during calculation.
  }
}
\details{
  This computes four standard summary statistics for a
  point pattern: the empty space function \eqn{F(r)},
  nearest neighbour distance distribution function \eqn{G(r)},
  van Lieshout-Baddeley function \eqn{J(r)}
  and Ripley's function \eqn{K(r)}.
  The real work is done by 
  \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}} and
  \code{\link{Kest}} respectively. Consult the help files for these functions
  for further information about the statistical interpretation
  of \eqn{F}, \eqn{G}, \eqn{J} and \eqn{K}.

  If \code{verb} is \code{TRUE}, then ``progress reports''
  (just indications of completion) are printed out when the
  calculations are finished for each of the four function types.

  The overall title of the array of four functions
  (for plotting by \code{\link{plot.fasp}})
  will be formed from the argument \code{dataname}.
  If this is not given, it defaults to the expression
  for \code{pp} given in the call to \code{allstats}.
}
\value{
  A function array (an object of class \code{"fasp"},
  see \code{\link{fasp.object}}). This can be plotted
  using \code{\link{plot.fasp}}.

  The function array has dimensions \eqn{2 \times 2}{2 * 2}
  with the entries at positions \code{[1,1]},
  \code{[1,2]}, \code{[2,1]} and \code{[2,2]} 
  representing \eqn{F(r)}, \eqn{G(r)}, \eqn{J(r)} and \eqn{K(r)} respectively.
  
  Each function entry \code{fns[[i]]} retains the format
  of the output of the relevant estimating routine
  \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}} or
  \code{\link{Kest}}. 

  The default formulae for plotting these functions are 
  \code{cbind(km,theo) ~ r} for F, G, and J, and
  \code{cbind(border,theo) ~ r} for K.}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
  }
\seealso{
         \code{\link{plot.fasp}},
         \code{\link{fasp.object}},
         \code{\link{Fest}},
         \code{\link{Gest}},
         \code{\link{Jest}},
         \code{\link{Kest}},
	 \code{\link{alltypes}}
}

\examples{
	data(swedishpines)
        a <- allstats(swedishpines,dataname="Swedish Pines")
        \dontrun{
        plot(a)
        plot(a, subset=list("r<=15","r<=15","r<=15","r<=50"))
        }
}
\keyword{spatial}
