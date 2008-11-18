#  tests/segments.R
#  $Revision: 1.4 $  $Date: 2008/11/16 06:16:26 $

require(spatstat)

# pointed out by Jeff Laake
W <- owin()
X <- psp(x0=.25,x1=.25,y0=0,y1=1,window=W)
X[W]

# test of distppll pointed out by Ang Qi Wei

p <- matrix(c(1.5, 0), 1, 2)
l <- matrix(c(0,0,1,0,1,0,2,0), 2, 4, byrow=T)
a <- distppll(p, l, mintype=2, method="interpreted")
b <- distppll(p, l, mintype=2, method="Fortran")
if(a$min.which != b$min.which)
  stop("conflict between Fortran and interpreted code in distppll")

