# check for various bugs related to factor conversions
require(spatstat)
# make a factor image
m <- factor(rep(letters[1:4], 4))
Z <- im(m, xcol=1:4, yrow=1:4)
# make a point pattern
set.seed(42)
X <- runifpoint(20, win=as.owin(Z))
# look up the image at the points of X
# (a) internal
ans1 <- lookup.im(Z, X$x, X$y)
stopifnot(is.factor(ans1))
# (b) user level
ans2 <- Z[X]
stopifnot(is.factor(ans2))
# (c) turn the image into a tessellation
#  and apply quadratcount
V <- tess(image = Z)
quadratcount(X, tess=V)
