#  tests/segments.R
#  $Revision: 1.5 $  $Date: 2010/04/20 11:54:26 $

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

# tests of pixellate.psp -> seg2pixL

ns <- 50
out <- numeric(ns)
for(i in 1:ns) {
  X <- psp(runif(1), runif(1), runif(1), runif(1), window=owin())
  len <- lengths.psp(X)
  dlen <- sum(pixellate(X)$v)
  out[i] <- if(len > 1e-7) dlen/len else 1
}
if(diff(range(out)) > 0.01) stop(paste(
       "pixellate.psp test 1: relative error [",
       paste(diff(range(out)), collapse=", "),
       "]"))

# Michael Sumner's test examples

set.seed(33)
n <- 2001
co <- cbind(runif(n), runif(n))
ow <- owin()
X <- psp(co[-n,1], co[-n,2], co[-1,1], co[-1,2], window=ow)
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 2:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

wts <- 1/(lengths.psp(X) * X$n)
s1 <- sum(pixellate(X, weights=wts))
if(abs(s1-1) > 0.01) {
  stop(paste("pixellate.psp test 3:",
             "sum(pixellate(X, weights))=", s1,
             " (should be 1)"))
}

X <- psp(0, 0, 0.01, 0.001, window=owin())
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 4:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

X <- psp(0, 0, 0.001, 0.001, window=owin())
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 5:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}




