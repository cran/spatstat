\name{nncorr}
\alias{nncorr}
\title{Nearest-Neighbour Correlation of Marked Point Pattern}
\description{
  Computes the nearest-neighbour correlation of a marked point pattern.
}
\usage{
     nncorr(X,
            f = function(m1, m2) { m1 * m2 },
            \dots,
            use = "all.obs", method = c("pearson", "kendall", "spearman"))
}
\arguments{
  \item{X}{
    The observed point pattern.
    An object of class \code{"ppp"}.
  }
  \item{f}{
    Function \eqn{f} used in the definition of the
    nearest neighbour correlation.
  }
  \item{\dots}{Extra arguments passed to \code{f}.}
  \item{use,method}{
    Arguments passed to the standard correlation function \code{\link{cor}}.
  }
}
\details{
  The nearest neighbour correlation \eqn{\bar n_f}{nbar}
  of a marked point process \eqn{X}
  is a number measuring the dependence between the mark of a typical point
  and the mark of its nearest neighbour.

  Three different values are computed: the unnormalised, normalised,
  and classical correlations.
  
  The \bold{unnormalised} nearest neighbour correlation (Stoyan and Stoyan,
  1994, section 14.7) is defined as
  \deqn{\bar n_f = E[f(M, M^\ast)]}{nbar[f] = E[f(M, M*)]}
  where \eqn{E[]} denotes mean value,
  \eqn{M} is the mark attached to a
  typical point of the point process, and \eqn{M^\ast}{M*} is the mark
  attached to its nearest neighbour (i.e. the nearest other point of the
  point process).
  
  Here \eqn{f} is any function
  \eqn{f(m_1,m_2)}{f(m1,m2)}
  with two arguments which are possible marks of the pattern,
  and which returns a nonnegative real value.
  Common choices of \eqn{f} are:
  for continuous real-valued marks,
  \deqn{f(m_1,m_2) = m_1 m_2}{f(m1,m2)= m1 * m2}
  for discrete marks (multitype point patterns),
  \deqn{f(m_1,m_2) = 1(m_1 = m_2)}{f(m1,m2)= (m1 == m2)}
  and for marks taking values in \eqn{[0,2\pi)}{[0,2 * pi)},
  \deqn{f(m_1,m_2) = \sin(m_1 - m_2)}{f(m1,m2) = sin(m1-m2).}
  For example, in the second case, the unnormalised nearest neighbour
  correlation \eqn{\bar n_f}{nbar[f]} equals the proportion of
  points in the pattern which have the same mark as their nearest
  neighbour.

  Note that \eqn{\bar n_f}{nbar[f]} is not a ``correlation''
  in the usual statistical sense. It can take values greater than 1.

  We can define a \bold{normalised} nearest neighbour correlation
  by 
  \deqn{\bar m_f = \frac{E[f(M,M^\ast)]}{E[f(M,M')]}}{mbar[f] = E[f(M,M*)]/E[f(M,M')]}
  where again \eqn{M} is the
  mark attached to a typical point, \eqn{M^\ast}{M*} is the mark
  attached to its nearest neighbour, and \eqn{M'} is an independent
  copy of \eqn{M} with the same distribution.
  This normalisation is also not a ``correlation''
  in the usual statistical sense, but is normalised so that 
  the value 1 suggests ``lack of correlation'':
  if the marks attached to the points of \code{X} are independent
  and identically distributed, then
  \eqn{\bar m_f = 1}{mbar[f] =  1}.
  The interpretation of values larger or smaller than 1 depends
  on the choice of function \eqn{f}.

  Finally if the marks of \code{X} are real numbers, we can also compute the
  \bold{classical} correlation, that is, the correlation coefficient
  of the two random variables \eqn{M} and \eqn{M^\ast}{M*}.
  The classical correlation has a value between \eqn{-1} and \eqn{1}.
  Values close to \eqn{-1} or \eqn{1} indicate strong dependence between
  the marks.

  This function computes the unnormalised and normalised
  nearest neighbour correlations, and the classical correlation
  if appropriate.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) and must be a marked point pattern.

  The argument \code{f}
  must be a function, accepting two arguments \code{m1}
  and \code{m2} which are vectors of equal length containing mark
  values (of the same type as the marks of \code{X}).
  It must return a vector of numeric
  values of the same length as \code{m1} and \code{m2}.
  The values must be non-negative.

  The arguments \code{use} and \code{method} control
  the calculation of the classical correlation using \code{\link{cor}},
  as explained in the help file for \code{\link{cor}}.

  Other arguments may be passed to \code{f} through the \code{...}
  argument.
  
  This algorithm assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape.
  Biases due to edge effects are
  treated using the \sQuote{border method} edge correction.
}
\value{
  Labelled vector of length 2 or 3 containing the unnormalised and normalised
  nearest neighbour correlations, and the classical correlation
  if appropriate.
}
\examples{
  data(finpines)
  nncorr(finpines)
  # heights of neighbouring trees are slightly negatively correlated

  data(amacrine)
  nncorr(amacrine, function(m1, m2) { m1 == m2})
  # neighbouring cells are usually of different type
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
