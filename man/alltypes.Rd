\name{alltypes}
\alias{alltypes}
\title{Calculate Statistic for All Types in a Multitype Point Pattern}
\description{
  Given a marked point pattern, this computes the estimates of
  a selected summary function (\eqn{F},\eqn{G}, \eqn{J} or \eqn{K})
  of the pattern, for all possible combinations of marks.
  It returns these functions in
  a list (an object of class \code{"fasp"}) amenable to plotting
  by \code{\link{plot.fasp}()}.
}

\usage{
  alltypes(pp, fun="K",dataname=NULL,verb=FALSE)
}

\arguments{
  \item{pp}{The observed point pattern, for which summary function
    estimates are required.  An object of class \code{"ppp"}.
    If the pattern is not marked, the resulting ``array'' has dimensions
    \eqn{1 \times 1}{1 x 1}.
  }

  \item{fun}{Character string indicating the summary function
    required.  Options are
    \code{"F"}, \code{"G"}, \code{"J"}, \code{"K"},
    \code{"Gcross"}, \code{"Jcross"}, \code{"Kcross"},
    \code{"Gdot"}, \code{"Jdot"}, \code{"Kdot"}.
  }

  \item{dataname}{Character string giving an optional (alternative)
    name to the point pattern, different from what is given
    in the call.  This name, if supplied, may be used by
    \code{\link{plot.fasp}()} in forming the title of the plot.
    If not supplied it defaults to the parsing of the argument
    supplied as \code{pp} in the call.
  }

  \item{verb}{Logical value, meaning ``verbose''.  If verb is
   true then terse ``progress reports'' (just the values of the
   mark indices) are printed out when the calculations for that
   combination of marks are completed. 
   }
}

\details{
  This routine is a convenient way to analyse the dependence between
  types in a multitype point pattern.
  Suppose that the points have possible types \eqn{1,2,\ldots,m}
  and let \eqn{X_i}{X[i]} denote the pattern of points of type \eqn{i} only.

  If \code{fun="F"} then this routine
  calculates, for each possible type \eqn{i},
  an estimate of the Empty Space Function \eqn{F_i(r)} of
  \eqn{X_i}{X[i]}. See \code{\link{Fest}} for explanation of the
  empty space function.

  If \code{fun} is
  \code{"Gcross"}, \code{"Jcross"} or \code{"Kcross"},
  the routine calculates, for each pair of types \eqn{(i,j)},
  an estimate of the ``\code{i}-to\code{j}'' cross-type function
  \eqn{G_{ij}(r)}{G[i,j](r)},
  \eqn{J_{ij}(r)}{J[i,j](r)} or
  \eqn{K_{ij}(r)}{K[i,j](r)} respectively describing the
  dependence between 
  \eqn{X_i}{X[i]} and \eqn{X_j}{X[j]}.
  See \code{\link{Gcross}}, \code{\link{Jcross}}
  or \code{\link{Kcross}} respectively for explanation of these functions.

  If \code{fun} is 
  \code{"Gdot"}, \code{"Jdot"} or \code{"Kdot"},
  the routine calculates, for each type \eqn{i},
  an estimate of the ``\code{i}-to-any'' dot-type function
  \eqn{G_{i\bullet}(r)}{G[i.](r)},
  \eqn{J_{i\bullet}(r)}{J[i.](r)} or
  \eqn{K_{i\bullet}(r)}{K[i.](r)} respectively describing the
  dependence between \eqn{X_i}{X[i]} and \eqn{X}{X}.
  See \code{\link{Gdot}}, \code{\link{Jdot}}
  or \code{\link{Kdot}} respectively for explanation of these functions.

  The letters \code{"G"}, \code{"J"} and \code{"K"}
  are interpreted as abbreviations for \code{"Gcross"}, \code{"Jcross"}
  and \code{"Kcross"} respectively, assumin the point pattern is
  marked. If the point pattern is unmarked, the appropriate
  function \code{\link{Fest}}, \code{\link{Jest}}
  or \code{\link{Kest}} is invoked instead.
  
}

\value{
  A function array (an object of class \code{"fasp"},
  see \code{\link{fasp.object}}). This can be plotted
  using \code{\link{plot.fasp}}.

  If \code{fun="F"},
  the function array has dimensions \eqn{m \times 1}{m * 1}
  where \eqn{m} is the number of different marks in the point pattern.
  The entry at position \code{[i,1]} in this array
  is the result of applying \code{\link{Fest}} to the
  points of type \code{i} only.

  If \code{fun} is \code{"Gdot"}, \code{"Jdot"}
  or \code{"Kdot"}, the function array
  again has dimensions \eqn{m \times 1}{m * 1}.
  The entry at position \code{[i,1]} in this array
  is the result of \code{Gdot(pp, i)} or \code{Jdot(pp, i)}
  or \code{Kdot(pp, i)} respectively.

  If \code{fun} is \code{"Gcross"}, \code{"Jcross"} or \code{"Kcross"},
  (or their abbreviations \code{"G"}, \code{"J"} or \code{"K"}),
  the function array has dimensions \eqn{m \times m}{m * m}.
  The \code{[i,j]} entry of the function array
  (for \eqn{i \neq j}{i != j}) is the
  result of applying the function \code{\link{Gcross}},
  \code{\link{Jcross}} or \code{\link{Kcross}} to
  the pair of types \code{(i,j)}. The diagonal
  \code{[i,i]} entry of the function array is the result of
  applying the univariate function \code{\link{Gest}},
  \code{\link{Jest}} or \code{\link{Kest}} to the
  points of type \code{i} only.
  
  Each function entry \code{fns[[i]]} retains the format
  of the output of the relevant estimating routine
  \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}},
  \code{\link{Kest}},  \code{Gcross}, \code{Jcross}, \code{Kcross},
  \code{Gdot}, \code{Jdot}, or \code{Kdot}.

  The default formulae for plotting these functions are 
  \code{cbind(km,theo) ~ r} for F, G, and J functions, and
  \code{cbind(trans,theo) ~ r} for K functions.
}
\note{
  Sizeable amounts of memory may be needed during the calculation.
}

\seealso{
  \code{\link{plot.fasp}},
  \code{\link{fasp.object}},
  \code{\link{allstats}},
  \code{\link{Fest}},
  \code{\link{Gest}},
  \code{\link{Jest}},
  \code{\link{Kest}},
  \code{\link{Gcross}},
  \code{\link{Jcross}},
  \code{\link{Kcross}},
  \code{\link{Gdot}},
  \code{\link{Jdot}},
  \code{\link{Kdot}}
}
\examples{
   # bramblecanes (3 marks).
   data(bramblecanes)
   \dontrun{
   X.F <- alltypes(bramblecanes,fun="F",verb=TRUE)
   plot(X.F) 
   X.G <- alltypes(bramblecanes,fun="G",verb=TRUE)
   X.J <- alltypes(bramblecanes,fun="J",verb=TRUE)
   X.K <- alltypes(bramblecanes,fun="K",verb=TRUE)
   X.Gd <- alltypes(bramblecanes,fun="Gdot",verb=TRUE)
   X.Jd <- alltypes(bramblecanes,fun="Jdot",verb=TRUE)
   X.Kd <- alltypes(bramblecanes,fun="Kdot",verb=TRUE)
   }
   \testonly{
      # smaller dataset
      bram <- bramblecanes[seq(1, bramblecanes$n, by=20), ]
      X.F <- alltypes(bram,fun="F",verb=TRUE)
      X.G <- alltypes(bram,fun="G",verb=TRUE)
      X.J <- alltypes(bram,fun="J",verb=TRUE)
      X.K <- alltypes(bram,fun="K",verb=TRUE)
      X.Gd <- alltypes(bram,fun="Gdot",verb=TRUE)
      X.Jd <- alltypes(bram,fun="Jdot",verb=TRUE)
      X.Kd <- alltypes(bram,fun="Kdot",verb=TRUE)
   }
   
   # Swedishpines (unmarked).
   data(swedishpines)
   \testonly{
	swedishpines <- swedishpines[1:25]
   }
   X.K <- alltypes(swedishpines,fun="K")
   X.F <- alltypes(swedishpines,fun="F")
   X.G <- alltypes(swedishpines,fun="G")
   X.J <- alltypes(swedishpines,fun="J")

   # simulated data
   \dontrun{
   pp <- runifpoint(350, owin(c(0,1),c(0,1)))
   pp$marks <- factor(c(rep(1,50),rep(2,100),rep(3,200)))
   X.F <- alltypes(pp,fun="F",verb=TRUE,dataname="Fake Data")
   X.G <- alltypes(pp,fun="G",verb=TRUE,dataname="Fake Data")
   X.J <- alltypes(pp,fun="J",verb=TRUE,dataname="Fake Data")
   X.K <- alltypes(pp,fun="K",verb=TRUE,dataname="Fake Data")
   }

   # A setting where you might REALLY want to use dataname:
   \dontrun{
   xxx <- alltypes(ppp(Melvin$x,Melvin$y,
                window=as.owin(c(5,20,15,50)),marks=clyde),
                fun="F",verb=TRUE,dataname="Melvin")
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
