# ppmmarkorder.R
# $Revision: 1.1 $  $Date: 2009/01/28 00:13:44 $
# Test that predict.ppm, plot.ppm and plot.fitin
# tolerate marks with levels that are not in alpha order
#
require(spatstat)
data(amacrine)
X <- amacrine
levels(marks(X)) <- c("ZZZ", "AAA")
fit <- ppm(X, ~marks, MultiStrauss(c("ZZZ","AAA"), matrix(0.06, 2, 2)))
aa <- predict(fit, type="trend")
bb <- predict(fit, type="cif")
plot(fit)
plot(fitin(fit))

