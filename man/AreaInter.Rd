\name{AreaInter}
\alias{AreaInter}
\title{The Area Interaction Point Process Model}
\description{
  Creates an instance of the Area Interaction point process model
  (marginal Widom-Rowlinson model) 
  which can then be fitted to point pattern data.
}
\usage{
  AreaInter(r)
}
\arguments{
  \item{r}{The radius of the discs in the area interaction process}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the area-interaction process with disc radius \eqn{r}.
}
\details{
  The (stationary) area interaction process
  (Baddeley and Van Lieshout, 1995) with disc radius \eqn{r},
  intensity parameter \eqn{\beta}{beta} and interaction parameter
  \eqn{\gamma}{gamma} is a point process with probability density
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \gamma^{-A(x)}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) gamma^(-A(x))
  }
  where 
  \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, and \eqn{A(x)} is the area of the region formed by
  the union of discs of radius \eqn{r} centred at the points
  \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]}.
  Here \eqn{\alpha}{alpha} is a normalising constant.

  The interaction parameter \eqn{\gamma}{gamma} can be any positive number.
  If \eqn{\gamma = 1}{gamma = 1} then the model reduces to a Poisson
  process. If \eqn{\gamma < 1}{gamma < 1} then the process is regular,
  while if \eqn{\gamma > 1}{gamma > 1} the process is clustered.
  The value \eqn{\gamma = 0}{eta = 0} corresponds to a hard core process
  with hard core radius \eqn{r} (interaction distance \eqn{2r}).

  In \pkg{spatstat} this model is reparametrised so that the interaction
  parameter is not affected by scale changes. In the new
  parametrisation, the density is 
  \deqn{
    f(x_1,\ldots,x_n) =
    \alpha \beta^{n(x)} \eta^{-B(x)}
  }{
    f(x_1,\ldots,x_n) =
    alpha . beta^n(x) eta^(-B(x))
  }
  where \eqn{\eta}{eta} is the new interaction parameter, and 
  \deqn{
    B(x) = \frac{A(x)}{\pi r^2}
  }{
    B(x) = A(x)/(pi * r^2)
  }
  is the normalised area (so that the discs have unit area).
  The value \eqn{\eta = 1}{eta = 1}
  again corresponds to a Poisson process.
  
  The \emph{nonstationary} area interaction process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location, rather than a constant beta. 
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the area interaction structure is
  yielded by the function \code{AreaInter()}. See the examples below.
 
  Note the only argument is the disc radius \code{r}.
  When \code{r} is fixed, the model becomes an exponential family.
  The canonical parameters \eqn{\log(\beta)}{log(beta)}
  and \eqn{\log(\eta)}{log(eta)}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{AreaInter()}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\section{Warnings}{
  The interaction distance of this process is equal to
  \code{2 * r}. Two discs of radius \code{r} overlap if their centres
  are closer than \code{2 * r} units apart.
}
\examples{
   AreaInter(r=0.1)
   # prints a sensible description of itself
   data(swedishpines)
   # Fit the stationary area interaction process to Swedish Pines data   
   ppm(swedishpines, ~1, AreaInter(r=2), rbord=4)
   data(cells) 
   # Fit the stationary area interaction process to `cells'
   ppm(cells, ~1, AreaInter(r=0.03), rbord=0.06)
   # Zero value of eta indicates hard core process.

   ppm(swedishpines, ~polynom(x/10,y/10,3), AreaInter(r=2), rbord=4)
   # fit a nonstationary area interaction with log-cubic polynomial trend
}
\references{
  Baddeley, A.J. and Van Lieshout, M.N.M. (1995).
  Area-interaction point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{47} (1995) 601--619.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{models}
