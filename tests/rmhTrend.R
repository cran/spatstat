#
# tests/rmhTrend.R
#
#  Problems with trend images (rmhmodel.ppm or rmhEngine)
#

require(spatstat)
set.seed(42)

# Bug folder 37 of 8 feb 2011
# rmhmodel.ppm -> predict.ppm
# + rmhResolveTypes -> is.subset.owin

data(demopat)
Z <- rescale(demopat, 7000)
X <- unmark(Z)
X1 <- split(Z)[[1]]
Int  <- density(X,dimyx=200)
Lint <- eval.im(log(npoints(X1)*Int/npoints(X)))
M    <- as.owin(Int)
MR   <- intersect.owin(M,expand.owin(M,0.5))
X1 <- X1[MR]
Fut  <- ppm(X1,~offset(Lint),covariates=list(Lint=Lint),
            inter=BadGey(r=c(0.03,0.05),sat=3))
Y   <- rmh(Fut,control=list(expand=M,nrep=1e3), verbose=FALSE)
