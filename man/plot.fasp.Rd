\name{plot.fasp}
\alias{plot.fasp}
\title{Plot a Function Array}
\description{
   Plots an array of summary functions, usually associated with a
   point pattern, stored in an object of class \code{"fasp"}.
   A method for \code{plot}.
}

\synopsis{
   \method{plot}{fasp}(x,formule=NULL,subset=NULL,lty=NULL,col=NULL,title=NULL,\dots, samex=TRUE)
}

\usage{
\method{plot}{fasp}(x)
\method{plot}{fasp}(x, formule=formule)
\method{plot}{fasp}(x, formule=formule, subset)
\method{plot}{fasp}(x, formule=formule, subset, lty, col, title, \dots, samex=TRUE)
}

\arguments{
  \item{x}{An object of class \code{"fasp"} representing a
    function array.
  }

  \item{formule}{A formula or list of formulae indicating what
    variables are to be plotted against what variable. Each formula is
    either an R language formula object, or a string that can be parsed
    as a formula. If
    \code{formule} is a list, its \eqn{k^{th}}{k-th} component
    should be applicable to the \eqn{(i,j)^{th}}{(i,j)-th}
    plot where \code{x$which[i,j]=k}.  If the formula is left
    as \code{NULL}, then \code{plot.fasp} attempts to use the component
    \code{default.formula} of \code{x}.  If that component is NULL
    as well, it gives up.

    In \bold{S-PLUS}, if this argument is given, then either it \bold{must}
    be supplied in the ``name=value'' format, \bold{or} the name
    of the plotting function \bold{must} be used explicitly.
    I.e. the call must use the function name \code{plot.fasp}
    and \bold{not} \code{plot}.  Otherwise the plotting function
    falls over.  See \bold{Warnings}.  There is no problem in R.
  }

  \item{subset}{A logical vector, or a vector of indices, or an
     expression or a character string, or a \bold{list} of such,
     indicating a subset of the data to be included in each plot.
     If \code{subset} is a list, its \eqn{k^{th}}{k-th} component
     should be applicable to the \eqn{(i,j)^{th}}{(i,j)-th} plot
     where \code{x$which[i,j]=k}.
  }

  \item{lty}{A vector of integers indicating line types to be
    used in forming the plot.
  }

  \item{col}{
    A vector of integers indicating colours for the different
    lines which form the plot.
  }

  \item{title}{
    Overall title for the plot.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.default}} to control other
    features of the individual plot panels.
  }
  \item{samex}{
    Logical flag indicating whether all individual plots should have the
    same x axis limits. This makes it easier to compare the plots.
    It can only be set to \code{FALSE} if you are using the default plot style
    (i.e. only when \code{formule} is missing).
  }
}

\details{
  An object of class \code{"fasp"} represents
  an array of summary functions, usually associated with a point
  pattern. See \code{\link{fasp.object}} for details.
  Such an object might be created, for example, by \code{\link{alltypes}}
  or \code{\link{allstats}}.

  The function \code{plot.fasp} is
  a method for \code{plot}.  It calls \code{\link{plot.fv}} to plot the
  individual panels.

  For information about the interpretation of the
  arguments \code{formule}, \code{subset}, \code{lty} and \code{col},
  see \code{\link{plot.fv}}.

  The argument \code{title}, if present, will determine the
  overall title of the plot. If it is absent, it defaults to \code{x$title}.
  Titles for the individual plot panels will be taken from
  \code{x$titles}.
}

\value{None.}

\section{Warnings}{
  (Each component of) the \code{subset} argument may be a
  logical vector (of the same length as the vectors of data which
  are extracted from \code{x}), or a vector of indices, or an
  \bold{expression} such as \code{expression(r<=0.2)}, or a text string,
  such as \code{"r<=0.2"}.

  Attempting a syntax such as \code{subset = r<=0.2} (without
  wrapping \code{r<=0.2} either in quote marks or in \code{expression()})
  will cause this function to fall over.

  Variables referred to in any formula must exist in the data frames
  stored in \code{x}.  What the names of these variables are will
  of course depend upon the nature of \code{x}.
}

\seealso{
  \code{\link{alltypes}},
  \code{\link{allstats}},
  \code{\link{plot.fv}},
  \code{\link{fasp.object}}
}

\examples{
   \dontrun{
   # Bramble Canes data.
   data(bramblecanes)

   X.G <- alltypes(bramblecanes,type="G",dataname="Bramblecanes",verb=TRUE)
   plot(X.G)
   plot(X.G,subset="r<=0.2")
   plot(X.G,formule=cbind(asin(sqrt(km)),
                    asin(sqrt(theo))) ~ asin(sqrt(theo)))
       plot(X.G,fo=cbind(km-theo,0)~r,"r<=0.2")

   # Swedish pines.
   data(swedishpines)
   X <- allstats(swedishpines,dataname="Swedish Pines")
   plot(X,subset=list("r<=20","r<=20","r<=20","r<=50"))

   # Simulated data.
   pp <- runifpoint(350, owin(c(0,1),c(0,1)))
   pp$marks <- factor(c(rep(1,50),rep(2,100),rep(3,200)))
   X.K <- alltypes(pp,type="K",verb=TRUE,dataname="Fake Data")
   plot(X.K,fo=cbind(border,theo)~theo,"theo<=0.75")
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
\keyword{hplot}
