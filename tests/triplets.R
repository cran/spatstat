# tests/triplets.R
# test code for triplet interaction
# $Revision: 1.2 $ $Date: 2011/11/07 07:48:55 $
require(spatstat)
fit <- ppm(redwood, ~1, Triplets(0.1))
fit
suffstat(fit)
# hard core
fithard <- ppm(cells, ~1, Triplets(0.05))
fithard
suffstat(fithard)
