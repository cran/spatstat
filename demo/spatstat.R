if(dev.cur() <= 1) {
  dd <- getOption("device")
  if(is.character(dd)) dd <- get(dd)
  dd()
}

oldpar <- par(ask = interactive() && dev.interactive(orNone=TRUE))
oldoptions <- options(warn=-1)

fanfare <- function(stuff) {
  plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
  text(0.5,0.5, stuff, cex=2.5)
}

par(mar=c(1,1,2,1)+0.1)

fanfare("Spatstat demonstration")

fanfare("I. Types of data")
data(swedishpines)
plot(swedishpines, main="Point pattern")

data(demopat)
plot(demopat, cols=c("green", "blue"), main="Multitype point pattern")

data(longleaf)
plot(longleaf, fg="blue", main="Marked point pattern")

a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
plot(a, main="Line segment pattern")
marks(a) <- sample(letters[1:4], 20, replace=TRUE)
plot(a, main="Multitype line segment pattern")
marks(a) <- runif(20)
plot(a, main="Marked line segment pattern")

plot(owin(), main="Rectangular window")
data(letterR)
plot(letterR, main="Polygonal window")
plot(as.mask(letterR), main="Binary mask window")

Z <- as.im(function(x,y){ sqrt((x - 1)^2 + (y-1)^2)}, square(2))
plot(Z, main="Pixel image")

X <- runifpoint(42)
plot(dirichlet(X), main="Tessellation")

enable3d <- ("scatterplot3d" %in% row.names(installed.packages()))
if(enable3d)
  plot(rpoispp3(100), main="Three-dimensional point pattern")

data(simplenet)
plot(simplenet, main="Linear network (linnet)")

X <- rpoislpp(20, simplenet)
plot(X, main="Point pattern on linear network (lpp)")

fanfare("II. Graphics")

plot(letterR, col="green", border="red", lwd=2, main="Polygonal window with colour fill")
plot(letterR, hatch=TRUE, spacing=0.15, angle=30, main="Polygonal window with line shading")

data(amacrine)
plot(amacrine, chars=c(1,16),
     main="plot(X, chars = c(1,16))")
plot(amacrine, cols=c("red","blue"), chars=16,
     main="plot(X, cols=c(\"red\", \"blue\"))")

opa <- par(mfrow=c(1,2))
plot(longleaf, markscale=0.03, main="markscale=0.03")
plot(longleaf, markscale=0.09, main="markscale=0.09")           
par(opa)
           
Z <- as.im(function(x,y) { r <- sqrt(x^2+y^2); r * exp(-r) },
           owin(c(-5,5),c(-5,5)))
plot(Z, main="pixel image: image plot")
plot(Z, main="pixel image: image plot (heat colours)", col=heat.colors(256))
contour(Z, main="pixel image: contour plot", axes=FALSE)
plot(Z, main="pixel image: image + contour plot")
contour(Z, add=TRUE)
persp(Z, colmap=terrain.colors(128), shade=0.3, phi=30,theta=100,
      main="pixel image: perspective plot")

ct <- colourmap(rainbow(20), breaks=seq(-1,1,length=21))
plot(ct, main="Colour map for real numbers")

ca <- colourmap(rainbow(8), inputs=letters[1:8])
plot(ca, main="Colour map for discrete values")

data(letterR)
W <- owin(c(1,5),c(0,4.5))
Lout <- scaletointerval(distmap(rebound.owin(letterR, W)))
Lin <- scaletointerval(distmap(complement.owin(letterR, W)))
L <- scaletointerval(eval.im(Lin-Lout))
D <- scaletointerval(density(runifpoint(30, W), adjust=0.3))
X <- scaletointerval(as.im(function(x,y){ x }, W=W))
plot(listof(L=L, D=D, X=X), main="Multiple images")
pairs(L, D, X, main="Multiple images: pairs plot")
plot(rgbim(D,X,L,maxColorValue=1), valuesAreColours=TRUE, main="Three images: RGB display")
plot(hsvim(D,L,X), valuesAreColours=TRUE, main="Three images: HSV display")

fanfare("III. Conversion between types")

data(chorley)
W <- as.owin(chorley)
plot(W, "window W")

plot(as.mask(W))
plot(as.mask(W, dimyx=1000))

plot(as.im(W, value=3))
plot(as.im(W, value=3, na.replace=0), ribbon=TRUE)

plot(as.im(function(x,y) {x^2 + y}, W=square(1)),
     main="as.im(function(x,y){x^2+y})")

V <- delaunay(runifpoint(12))
plot(V, main="Tessellation V")
plot(as.im(V, dimyx=256), main="as.im(V)")
plot(as.owin(V))

X <- swedishpines
plot(X, "point pattern X")

plot(as.im(X), col=c("white","red"), ribbon=FALSE, xlab="", ylab="")
plot(as.owin(X), add=TRUE)

fanfare("IV. Subsetting and splitting data")

plot(X, "point pattern X")
subset <- 1:20
plot(X[subset], main="subset operation: X[subset]")
subwindow <- owin(poly=list(x=c(0,96,96,40,40),y=c(0,0,100,100,50)))
plot(X[subwindow], main="subset operation: X[subwindow]")

data(lansing)
plot(lansing, "Lansing Woods data")
plot(split(lansing), main="split operation: split(X)")

data(longleaf)
plot(longleaf, main="Longleaf Pines data")
plot(cut(longleaf, breaks=3),
     main=c("cut operation", "cut(longleaf, breaks=3)"))

Z <- dirichlet(runifpoint(16))
X <- runifpoint(100)

plot(Z, main="points cut by tessellation")
plot(cut(X, Z), add=TRUE)

plot(split(X, Z), main="points split by tessellation")

W <- square(1)
X <- as.im(function(x,y){sqrt(x^2+y^2)}, W)
Y <- dirichlet(runifpoint(12, W))
plot(split(X,Y), main="image split by tessellation")

fanfare("V. Exploratory data analysis")

par(mar=c(3,3,3,2)+0.1)

plot(swedishpines, main="Quadrat counts", pch="+")
tab <- quadratcount(swedishpines, 4)
plot(tab, add=TRUE, lty=2, cex=2, col="blue")

par(mar=c(5,3,3,2)+0.1)

plot(swedishpines, main="", pch="+")
title(main=expression(chi^2 * " test"), cex.main=2)
tes <- quadrat.test(swedishpines, 3)
tes
plot(tes, add=TRUE, col="red", cex=1.5, lty=2, lwd=3)
title(sub=paste("p-value =", signif(tes$p.value,3)), cex.sub=1.4)

par(mar=c(4,4,3,2)+0.1)

data(nztrees)
tesk <- kstest(nztrees, "x")
tesk
plot(tesk)

data(murchison)
mur <- lapply(murchison, rescale, s=1000)
X <- mur$gold
D <- distfun(mur$faults)
rh <- rhohat(X,D, dimyx=256)
plot(rh,
     main="Smoothed rate estimate",
     xlab="Distance to nearest fault (km)",
     legend=FALSE)
plot(predict(rh), main="predict(rhohat(X,D))")

data(cells)
Z <- density.ppp(cells, 0.07)
plot(Z, main="Kernel smoothed intensity of point pattern")
plot(cells, add=TRUE)

data(shapley)
X <- unique(unmark(shapley))
plot(X, "Shapley galaxy concentration", pch=".")
plot(nnclean(X, k=17), main="Byers-Raftery nearest neighbour cleaning",
     chars=c(".", "+"), cols=1:2)
Y <- sharpen(X, sigma=0.5, edgecorrect=TRUE)
plot(Y, main="Choi-Hall data sharpening", pch=".")

D <- density(a, sigma=0.05)
plot(D, main="Kernel smoothed intensity of line segment pattern")
plot(a, add=TRUE)

X <- runifpoint(42)
plot(dirichlet(X))
plot(X, add=TRUE)

plot(delaunay(X))
plot(X, add=TRUE)

data(longleaf)
parsave <- par(mfrow=c(1,2))
plot(longleaf, main="Longleaf Pines data")
plot(smooth.ppp(longleaf, 10), main="Spatial smoothing of marks")
par(parsave)

data(cells)
fryplot(cells, main=c("Fry plot","cells data"), pch="+")
data(longleaf)
miplot(longleaf, main="Morishita Index plot", pch=16, col="blue")

plot(swedishpines, main="Swedish Pines data")
K <- Kest(swedishpines)
plot(K, main="K function for Swedish Pines", legendmath=TRUE)

en <- envelope(swedishpines, fun=Kest, nsim=10, correction="translate")
plot(en, main="Envelopes of K function based on CSR", shade=c("hi", "lo"))

pc <- pcf(swedishpines)
plot(pc, main="Pair correlation function")

plot(swedishpines, main="nearest neighbours")
m <- nnwhich(swedishpines)
b <- swedishpines[m]
arrows(swedishpines$x, swedishpines$y, b$x, b$y,
       angle=12, length=0.1, col="red")

plot(swedishpines %mark% (nndist(swedishpines)/2), markscale=1, main="Stienen diagram")

plot(Gest(swedishpines),
     main=c("Nearest neighbour distance function G", "Gest(swedishpines)"),
     legendmath=TRUE)

Z <- distmap(swedishpines, dimyx=512)
plot(swedishpines$window, main="Distance map")
plot(Z, add=TRUE)
points(swedishpines)

plot(Fest(swedishpines),
     main=c("Empty space function F", "Fest(swedishpines)"),
     legendmath=TRUE)

W <- rebound.owin(letterR, square(5))
plot(distmap(W), main="Distance map")
plot(W, add=TRUE)

a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
contour(distmap(a), main="Distance map")
plot(a, add=TRUE,col="red")

plot(Jest(swedishpines), main=c("J-function", "J(r)=(1-G(r))/(1-F(r))"))
     
plot(allstats(swedishpines))

data(residualspaper)
Fig4b <- residualspaper$Fig4b

plot(Fig4b, main="Inhomogeneous point pattern")
plot(Kinhom(Fig4b), main="Inhomogeneous K-function")
plot(pcfinhom(Fig4b, stoyan=0.1), main="Inhomogeneous pair correlation")

data(bronzefilter)
X <- unmark(bronzefilter)
plot(X, "Bronze filter data")
lam <- predict(ppm(X, ~x))
plot(Kscaled(X, lam), xlim=c(0, 1.5), main="Locally-scaled K function")

data(bramblecanes)
plot(bramblecanes)
bramblecanes <- rescale(bramblecanes, 1/9)
plot(alltypes(bramblecanes, "K"), mar.panel=c(4,5,2,2)+0.1)

data(amacrine)
amacrine <- rescale(amacrine, 1/662)
plot(alltypes(amacrine, Lcross, envelope=TRUE, nsim=9), . - r ~ r, ylim=c(-25, 5))

data(ponderosa)
ponderosa.extra$plotit(main="Ponderosa Pines")

L <- localL(ponderosa)
pL <- plot(L, lty=1, col=1, legend=FALSE,
           main="neighbourhood density functions for Ponderosa Pines")

parsave <- par(mfrow=c(1,2))
ponderosa.extra$plotit()
par(pty="s")
plot(L, iso007 ~ r, main="point B")

ponderosa.extra$plotit()
L12 <- localL(ponderosa, rvalue=12)
P12 <- ponderosa %mark% L12
Z12 <- smooth.ppp(P12, sigma=5, dimyx=128)
plot(Z12, col=topo.colors(128), main="smoothed neighbourhood density")
contour(Z12, add=TRUE)
points(ponderosa, pch=16, cex=0.5)

plot(amacrine, main="Amacrine cells data")
par(pty="s")
mkc <- markcorr(amacrine, 
                correction="translate", method="density",
                kernel="epanechnikov")
plot(mkc, main="Mark correlation function", legend=FALSE)
par(parsave)

plot(alltypes(amacrine, markconnect), 
     title="Mark connection functions for amacrine cells")

parsave <- par(mfrow=c(1,2))
data(spruces)
plot(spruces, cex.main=0.75)
par(pty="s")
plot(markcorr(spruces), main="Mark correlation", legendpos="bottomright")

plot(spruces, cex.main=0.75)
plot(markvario(spruces), main="Mark variogram", legendpos="topright")
par(parsave)

plot(as.listof(list("Emark(spruces)"=Emark(spruces),
                    "Vmark(spruces)"=Vmark(spruces))),
     main="Independence diagnostics", ylim.covers=0, legendpos="bottom")

if(enable3d) {
  par3 <- par(mfrow=c(1,2))
  X <- rpoispp3(100)
  plot(X, main="3D point pattern X")
  plot(K3est(X), main="K-function in 3D")
  plot(X, main="3D point pattern X")
  plot(G3est(X), main="G-function in 3D", legendpos="bottomright")
  par(par3)
}

par2 <- par(mfrow=c(1,3))
data(chicago)
X <- unmark(chicago)
plot(as.linnet(X), main="Chicago Street Crimes",col="green")
plot(as.ppp(X), add=TRUE, col="red")
plot(linearK(X, correction="none"), main="Network K-function")
plot(linearK(X, correction="Ang"), main="Corrected K-function")
par(par2)

fanfare("VI. Model-fitting")

data(japanesepines)
plot(japanesepines)
fit <- ppm(japanesepines, ~1)
print(fit)
fit <- ppm(japanesepines, ~polynom(x,y,2))
print(fit)
plot(fit, how="image", se=FALSE, main=c("Inhomogeneous Poisson model",
                               "fit by maximum likelihood",
                               "Fitted intensity"))
plot(fit, how="image", trend=FALSE,
     main="Standard error of fitted intensity")

data(redwood)
parsave <- par(mfrow=c(1,2))
plot(redwood)
fitT <- kppm(redwood, ~1, clusters="Thomas")
oop <- par(pty="s")
plot(fitT, main=c("Thomas model","fit by minimum contrast"))

plot(redwood)
plot(simulate(fitT)[[1]], main="simulation from fitted Thomas model")

plot(swedishpines)
fit <- ppm(swedishpines, ~1, Strauss(r=7))
print(fit)
plot(fit, how="image", main=c("Strauss model",
                               "fit by maximum pseudolikelihood",
                               "Conditional intensity plot"))
plot(Kcom(fit), cbind(iso, icom, pois) ~ r,
     legend=FALSE, main="model compensators")
legend("topleft", legend=c("empirical K function",
                    "Strauss model compensator of K",
                    "Poisson theoretical K"), lty=1:3, col=1:3, inset=0.05)

plot(swedishpines)
fit <- ppm(swedishpines, ~1, PairPiece(c(3,5,7,9,11,13)))
plot(fitin(fit), legend=FALSE,
     main=c("Pairwise interaction model",
            "fit by maximum pseudolikelihood"))
par(parsave)

Xsim <- rmh(model=fit,
            start=list(n.start=80),
            control=list(nrep=100))
plot(Xsim, main="Simulation from fitted Strauss model")


data(demopat)
demopat <- rescale(demopat, 8)
unitname(demopat) <- c("mile", "miles")
demopat

plot(demopat, cols=c("red", "blue"))
fit <- ppm(demopat, ~marks + polynom(x,y,2), Poisson())
plot(fit, trend=TRUE, se=TRUE)

fanfare("VII. Simulation")

data(letterR)
plot(letterR, main="Poisson random points")
lambda <- 10/area.owin(letterR)
points(rpoispp(lambda, win=letterR))
points(rpoispp(9 * lambda, win=letterR))
points(rpoispp(90 * lambda, win=letterR))
plot(rpoispp(100))
plot(rpoispp(function(x,y){1000 * exp(-3*x)}, 1000))

plot(rMaternII(200, 0.05))
plot(rSSI(0.05, 200))
plot(rThomas(10, 0.2, 5))
plot(rMatClust(10, 0.05, 4))
plot(rGaussPoisson(30, 0.05, 0.5))
param <- c(0, variance=0.2, nugget=0, scale=.1)
mu <- 4
plot(rLGCP("exp", mu, param))
X <- rLGCP("exp", mu, param)
plot(attr(X, "Lambda"), main="log-Gaussian Cox process")
plot(X, add=TRUE, pch=16)

plot(redwood, main="random thinning - rthin()")
points(rthin(redwood, 0.5), col="green", cex=1.4)

plot(rcell(nx=15))

plot(rsyst(nx=5))
abline(h=(1:4)/5, lty=2)
abline(v=(1:4)/5, lty=2)

plot(rstrat(nx=5))
abline(h=(1:4)/5, lty=2)
abline(v=(1:4)/5, lty=2)

X <- rsyst(nx=10)
plot(rjitter(X, 0.02))

Xg <- rmh(list(cif="geyer", par=list(beta=1.25, gamma=1.6, r=0.2, sat=4.5),
               w=c(0,10,0,10)),
          control=list(nrep=1e4), start=list(n.start=200))
plot(Xg, main=paste("Geyer saturation process\n",
                    "rmh() with cif=\"geyer\""))

L <- as.psp(matrix(runif(20), 5, 4), window=square(1))
plot(L, main="runifpointOnLines(30, L)")
plot(runifpointOnLines(30, L), add=TRUE, pch="+")

plot(L, main="rpoisppOnLines(3, L)")
plot(rpoisppOnLines(3, L), add=TRUE, pch="+")

data(simplenet)
plot(runiflpp(20, simplenet))
plot(rpoislpp(5, simplenet))

plot(rpoisline(10))

plot(rlinegrid(30, 0.1))

spatstat.options(npixel=256)
X <- dirichlet(runifpoint(30))
plot(rMosaicSet(X, 0.4), col="green", border=NA)
plot(X, add=TRUE)
plot(rMosaicField(X, runif))
plot(rMosaicSet(rpoislinetess(3), 0.5), col="green", border=NA, main="Switzer's random set")
spatstat.options(npixel=100)

fanfare("VIII. Geometry")

data(letterR)
A <- letterR

B <- shift(letterR, c(0.2,0.1))
plot(bounding.box(A,B), main="shift", type="n")
plot(A, add=TRUE)
plot(B, add=TRUE, border="red")

B <- rotate(letterR, 0.2)
plot(bounding.box(A,B), main="rotate", type="n")
plot(A, add=TRUE)
plot(B, add=TRUE, border="red")

mat <- matrix(c(1.1, 0, 0.3, 1), 2, 2)
B <- affine(letterR, mat=mat, vec=c(0.2,-0.1))
plot(bounding.box(A,B), main="affine", type="n")
plot(A, add=TRUE)
plot(B, add=TRUE, border="red")

par1x2 <- par(mfrow=c(1,2))
L <- rpoisline(10, owin(c(1.5,4.5),c(0.2,3.6)))
plot(L, main="Line segment pattern")
plot(L$window, main="L[window]", type="n")
plot(L[letterR], add=TRUE)
plot(letterR, add=TRUE, border="red")
par(par1x2)

a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
plot(a, main="Self-crossing points")
plot(selfcrossing.psp(a), add=TRUE, col="red")

a <- as.psp(matrix(runif(20), 5, 4), window=square(1))
b <- rstrat(square(1), 5)
plot(a, lwd=3, col="green", main="project points to segments")
plot(b, add=TRUE, col="red", pch=16)
v <- project2segment(b, a)
Xproj <- v$Xproj
plot(Xproj, add=TRUE, pch=16)
arrows(b$x, b$y, Xproj$x, Xproj$y, angle=10, length=0.15, col="red")

plot(a, main="pointsOnLines(L)")
plot(pointsOnLines(a, np=100), add=TRUE, pch="+")

parry <- par(mfrow=c(1,3))
X <- tess(xgrid=seq(2, 4, length=10), ygrid=seq(0, 3.5, length=8))
plot(X)
data(letterR)
plot(letterR)
plot(intersect.tess(X, letterR))

X <- dirichlet(runifpoint(10))
plot(X)
L <- infline(0.3,0.5)
plot(owin(), main="L")
plot(L, col="red", lwd=2)
plot(chop.tess(X,L))
par(parry)

data(chorley)
W <- chorley$window
plot(W, main="simplify.owin")
WS <- simplify.owin(W, 2)
plot(WS, add=TRUE, border="green")

nopa <- par(mfrow=c(2,2))
data(letterR)
Rbox <- grow.rectangle(as.rectangle(letterR), 0.3)

v <- erode.owin(letterR, 0.25)
plot(Rbox, type="n", main="erode.owin", cex.main=0.75)
plot(letterR, add=TRUE, col="red", cex.main=0.75)
plot(v, add=TRUE, col="blue")

v <- dilate.owin(letterR, 0.25)
plot(Rbox, type="n", main="dilate.owin", cex.main=0.75)
plot(v, add=TRUE, col="blue")
plot(letterR, add=TRUE, col="red")

v <- closing.owin(letterR, 0.3)
plot(Rbox, type="n", main="closing.owin", cex.main=0.75)
plot(v, add=TRUE, col="blue")
plot(letterR, add=TRUE, col="red")

v <- opening.owin(letterR, 0.3)
plot(Rbox, type="n", main="opening.owin", cex.main=0.75)
plot(letterR, add=TRUE, col="red")
plot(v, add=TRUE, col="blue")
par(nopa)

fanfare("IX. Operations on pixel images")

plot(Z, main="An image Z")
plot(levelset(Z, 4))
plot(cut(Z, 5))
plot(eval.im(sqrt(Z) - 3))
plot(solutionset(abs(Z - 6) <= 1))

data(cells)
d <- distmap(cells, dimyx=256)
W <- levelset(d, 0.06)
nopa <- par(mfrow=c(1,2))
plot(W)
plot(connected(W))
par(nopa)

Z <- as.im(function(x,y) { 4 * x^2 + 3 * y }, letterR)
plot(Z)
plot(letterR, add=TRUE)

plot(blur(Z, 0.3, bleed=TRUE))
plot(letterR, add=TRUE)

plot(blur(Z, 0.3, bleed=FALSE))
plot(letterR, add=TRUE)
          
plot(blur(Z, 0.3, bleed=FALSE))
plot(letterR, add=TRUE)
          

fanfare("X. Programming tools")

showoffK <- function(Y, current, ..., fullpicture,rad) { 
	plot(fullpicture,
             main=c("Animation using `applynbd'", "explaining the K function"))
	points(Y, cex=2)
        u <- current
	points(u[1],u[2],pch="+",cex=3)
	theta <- seq(0,2*pi,length=100)
	polygon(u[1]+ rad * cos(theta),u[2]+rad*sin(theta))
	text(u[1]+rad/3,u[2]+rad/2,Y$n,cex=3)
        if(runif(1) < 0.2) Sys.sleep(runif(1, max=0.4))
	return(Y$n)
}
par(ask=FALSE)
applynbd(redwood, R=0.2, showoffK, fullpicture=redwood, rad=0.2, exclude=TRUE)

par(oldpar)
options(oldoptions)

