# temporary test file for localpcfmatrix

require(spatstat)
data(redwood)
a <- localpcfmatrix(redwood)
a
plot(a)
a[, 3:5]
