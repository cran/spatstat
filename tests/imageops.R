#
#  tests/imageops.R
#
#   $Revision: 1.4 $   $Date: 2011/10/04 05:47:35 $
#

require(spatstat)
A <- as.im(owin())
B <- as.im(owin(c(1.1, 1.9), c(0,1)))
Z <- imcov(A, B)
stopifnot(abs(max(Z) - 0.8) < 0.1)



