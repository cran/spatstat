# ppmBadData.R
# $Revision: 1.3 $ $Date: 2011/01/24 07:00:11 $

# Testing robustness of ppm and support functions
# when data are rubbish

require(spatstat)

# ---------------------------------------------------
# from Rolf: very large proportion of data is NA

SEED <- 42

K <- 101
A <- 500
X <- seq(0, A, length=K)
G <- expand.grid(x=X, y=X)
FOO <- function(x,y) { sin(x)^2 + cos(y)^2 }
M1 <- im(matrix(FOO(G$x, G$y), K, K), xcol=X, yrow=X)
M <- im(matrix(FOO(G$x, G$y), K, K))
BAR <- function(x) { exp(-6.618913 + 5.855337 * x - 8.432483 * x^2) }
V <- im(BAR(M$v), xcol=X, yrow=X)
# V <- eval.im(exp(-6.618913 + 5.855337 * M - 8.432483 * M^2))
set.seed(SEED)
Y <- rpoispp(V)
fY <- ppm(Y, ~cv + I(cv^2), covariates=list(cv=M), correction="translate")
diagnose.ppm(fY)
lurking(fY, covariate=as.im(function(x,y){x}, square(A)), type="raw")

# --------------------------------------------------------
# from Andrew Bevan: numerical overflow, ill-conditioned Fisher information

nongranite<- owin(poly = list(x = c(0, 8500, 7000, 6400, 6400, 6700, 7000, 7200, 7300, 8000, 8100, 8800, 9500, 10000, 10000, 0), y = c(0, 0, 2000, 3800, 4000, 5000, 6500, 7400, 7500, 8000, 8100, 9000, 9500, 9600, 10000, 10000)))

#Trend on raster grid
rain <- as.im(X=function(x,y) { x^2 + y^2 }, W=nongranite, dimyx=100)

#Generate a point pattern via a Lennard-Jones process
set.seed(SEED)
mod4<- rmhmodel(cif="lennard",
                par=list(beta=1, sigma=250, epsilon=2.2),
                trend=rain, w=nongranite)
ljtr<- rmh(mod4, start=list(n.start=80), control=list(p=1, nrep=1e5))

#Fit a point process model to the pattern with rain as a covariate
# NOTE INCORRECT TREND FORMULA
ljtrmod <- ppm(ljtr, trend= ~ Z, interaction=NULL, covariates=list(Z=rain))
ss <- summary(ljtrmod)

