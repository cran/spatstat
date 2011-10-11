#
# tests/rmhExpand.R
#
# test decisions about expansion of simulation window
#
#  $Revision: 1.1 $  $Date: 2011/10/07 04:08:59 $
#

require(spatstat)
data(cells)
fit <- ppm(cells, ~x)

# check rmhmodel.ppm
mod <- rmhmodel(fit)
wsim <- as.rectangle(mod$trend)
if(!identical(wsim, as.owin(cells)))
  stop("Expansion occurred improperly in rmhmodel.ppm")


