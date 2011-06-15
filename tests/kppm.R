#
# tests/kppm.R
#
# $Revision: 1.3 $ $Date: 2011/06/08 11:50:10 $
#
# Test functionality of kppm that depends on RandomFields
#

require(spatstat)
data(redwood)

fit0 <- kppm(redwood, ~1, "LGCP")
simulate(fit0)

fit <- kppm(redwood, ~x, "LGCP", covmodel=list(model="matern", nu=0.3))
simulate(fit)



