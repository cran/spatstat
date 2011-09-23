#       
#        tests/hobjects.R
#
#   Validity of methods for ppm(... method="ho")
#

require(spatstat)

data(cells)

fit  <- ppm(cells, ~1,         Strauss(0.1), method="ho", nsim=10)
fitx <- ppm(cells, ~offset(x), Strauss(0.1), method="ho", nsim=10)

a  <- AIC(fit)
ax <- AIC(fitx)

f  <- fitted(fit)
fx <- fitted(fitx)

p  <- predict(fit)
px <- predict(fitx)

