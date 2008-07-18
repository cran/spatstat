# test for step() operation
#
require(spatstat)
data(nztrees)
Z <- as.im(function(x,y){ x^3 - y^2 }, nztrees$window)
fitP <- ppm(nztrees, ~x+y+Z, covariates=list(Z=Z))
step(fitP)
fitS <- update(fitP, Strauss(7))
step(fitS)
data(amacrine)
fitM <- ppm(amacrine, ~ marks*(x+y),
            MultiStrauss(types=levels(marks(amacrine)), radii=matrix(0.04, 2, 2)))
step(fitM)

