#  tests/segments.R
#  $Revision: 1.3 $  $Date: 2008/04/16 13:26:25 $

require(spatstat)

# pointed out by Jeff Laake
W <- owin()
X <- psp(x0=.25,x1=.25,y0=0,y1=1,window=W)
X[W]

