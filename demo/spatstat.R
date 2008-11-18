if(dev.cur() <= 1) get(getOption("device"))()

oldpar <- par(ask = interactive() && dev.interactive(orNone=TRUE))
oldoptions <- options(warn=-1)

fanfare <- function(stuff) {
  plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
  text(0.5,0.5, stuff, cex=2.5)
}
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

plot(owin(), main="Rectangular window")
data(letterR)
plot(letterR, main="Polygonal window")
plot(as.mask(letterR), main="Binary mask window")

plot(letterR, col="green", border="red", lwd=2, main="Polygonal window with colour fill")
plot(letterR, hatch=TRUE, spacing=0.15, angle=30, main="Polygonal window with line shading")

Z <- as.im(function(x,y){ sqrt((x - 1)^2 + (y-1)^2)}, square(2))
plot(Z, main="Pixel image")
plot(Z, main="Pixel image", col=heat.colors(256))

X <- runifpoint(42)
plot(dirichlet(X), main="Tessellation")

fanfare("II. Basic operations")
X <- swedishpines
subset <- 1:20
plot(X[subset], main="subset operation: X[subset]")
subwindow <- owin(poly=list(x=c(0,96,96,40,40),y=c(0,0,100,100,50)))
plot(X[subwindow], main="subset operation: X[subwindow]")

L <- rpoisline(10, owin(c(1.5,4.5),c(0.2,3.6)))
S <- L[letterR]
plot(L, main="subset operation: L[subwindow]")
plot(S, add=TRUE, col="red")

data(lansing)
plot(lansing, "Lansing Woods data")
plot(split(lansing), main="split operation: split(X)")

data(longleaf)
plot(longleaf, main="Longleaf Pines data")
plot(cut(longleaf, breaks=3),
     main=c("cut by marks", "cut(longleaf, breaks=3)"))

X <- runifpoint(100)
Z <- dirichlet(runifpoint(16))
plot(Z, main="cut by tessellation")
plot(cut(X, Z), add=TRUE)

plot(split(X, Z), main="split by tessellation")

W <- square(1)
X <- as.im(function(x,y){sqrt(x^2+y^2)}, W)
Y <- dirichlet(runifpoint(12, W))
plot(split(X,Y), main="image split by tessellation")

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

fanfare("III. Exploratory data analysis")

plot(swedishpines, main="Quadrat counts", pch="+")
tab <- quadratcount(swedishpines, 4)
plot(tab, add=TRUE, lty=2, cex=2, col="blue")

plot(swedishpines, main="", pch="+")
title(main=expression(chi^2 * " test"), cex.main=2)
tes <- quadrat.test(swedishpines, 3)
tes
plot(tes, add=TRUE, col="red", cex=1.5, lty=2, lwd=3)
title(sub=paste("p-value =", signif(tes$p.value,3)), cex.sub=1.4)

tesk <- ks.test.ppm(ppm(swedishpines), function(x,y){x})
tesk
plot(tesk)

data(cells)
Z <- density.ppp(cells, 0.07)
plot(Z, main="Kernel smoothed intensity of point pattern")
plot(cells, add=TRUE)

data(bei)
ZA <- adaptive.density(bei, 0.01, nrep=5)
plot(ZA, main="Adaptive intensity of point pattern",
     col=grey(seq(1,0,length=256)))
plot(bei, add=TRUE, pch=".")

D <- density(a, sigma=0.05)
plot(D, main="Kernel smoothed intensity of line segment pattern")
plot(a, add=TRUE)

data(longleaf)
parsave <- par(mfrow=c(1,2))
plot(longleaf, main="Longleaf Pines data")
plot(smooth.ppp(longleaf, 10), main="Spatial smoothing of marks")
par(parsave)

plot(swedishpines, main="Swedish Pines data")
K <- Kest(swedishpines)
plot(K, main="K function for Swedish Pines")

en <- envelope(swedishpines, fun=Kest, nsim=10, correction="translate")
plot(en, main="Envelopes of K function based on CSR")

pc <- pcf(swedishpines)
plot(pc, main="Pair correlation function")


plot(swedishpines, main="nearest neighbours")
m <- nnwhich(swedishpines)
b <- swedishpines[m]
arrows(swedishpines$x, swedishpines$y, b$x, b$y,
       angle=12, length=0.1, col="red")

plot(swedishpines %mark% (nndist(swedishpines)/2), markscale=1, main="Stienen diagram")

Z <- distmap(swedishpines, dimyx=512)
plot(swedishpines$window, main="Distance map")
plot(Z, add=TRUE)
points(swedishpines)

W <- rebound.owin(letterR, square(5))
plot(distmap(W), main="Distance map")
plot(W, add=TRUE)

persp(Z, colmap=terrain.colors(128), shade=0.3, phi=30,theta=100,
      main="perspective plot of pixel image")

a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
contour(distmap(a), main="Distance map")
plot(a, add=TRUE,col="red")

plot(allstats(swedishpines))

data(bramblecanes)
plot(bramblecanes)
bramblecanes <- rescale(bramblecanes, 1/9)
plot(alltypes(bramblecanes, "K"), ylab="K(r)", mar.panel=c(4,4,2,2)+0.1)

data(amacrine)
plot(alltypes(amacrine, Lcross, envelope=TRUE, nsim=9), ylab="L(r)")

data(ponderosa)
ponderosa.extra$plotit(main="Ponderosa Pines")

L <- localL(ponderosa)
plot(L, lty=1, col=1,
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

data(amacrine)
plot(amacrine, main="Amacrine cells data")
par(pty="s")
mkc <- markcorr(amacrine, function(m1,m2) {m1==m2},
                correction="translate", method="density",
                kernel="epanechnikov")
plot(mkc, main="Mark correlation function")

par(parsave)

X <- runifpoint(42)
plot(dirichlet(X))
plot(X, add=TRUE)

plot(delaunay(X))
plot(X, add=TRUE)

fanfare("IV. Model-fitting")

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
par(parsave)
plot(simulate(fitT)[[1]], main="simulation from fitted Thomas model")

plot(swedishpines)
fit <- ppm(swedishpines, ~1, Strauss(r=7))
print(fit)
plot(fit, how="image", main=c("Strauss model",
                               "fit by maximum pseudolikelihood",
                               "Conditional intensity plot"))

plot(swedishpines)
fit <- ppm(swedishpines, ~1, PairPiece(c(3,5,7,9,11,13)))
plot(fitin(fit),
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

fanfare("V. Simulation")

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

Xg <- rmh(list(cif="geyer", par=c(beta=1.25, gamma=1.6, r=0.2, sat=4.5),
               w=c(0,10,0,10)),
          control=list(nrep=1e4), start=list(n.start=200))
plot(Xg, main=paste("Geyer saturation process\n",
                    "rmh() with cif=\"geyer\""))

plot(rpoisline(10))

plot(rlinegrid(30, 0.1))

fanfare("VI. Programming tools")

nopa <- par(mfrow=c(2,2))
data(letterR)
Rbox <- as.rectangle(letterR)
Rmask <- as.mask(letterR, dimyx=256)

v <- erode.owin(Rmask, 0.25)
plot(Rbox, type="n", main="erode.owin")
plot(v, add=TRUE)
plot(letterR, add=TRUE)

v <- dilate.owin(w, 0.3)
plot(as.rectangle(v), type="n", main="dilate.owin")
plot(v, add=TRUE)
plot(letterR, add=TRUE)

v <- closing.owin(w, 0.25)
plot(Rbox, type="n", main="closing.owin")
plot(v, add=TRUE)
plot(letterR, add=TRUE)

v <- opening.owin(w, 0.3)
plot(Rbox, type="n", main="opening.owin")
plot(v, add=TRUE)
plot(letterR, add=TRUE)
par(nopa)

plot(Z, main="An image Z")
plot(levelset(Z, 4))
plot(cut(Z, 5))
plot(eval.im(sqrt(Z) - 3))
plot(solutionset(abs(Z - 6) <= 1))

Z <- as.im(function(x,y) { 4 * x^2 + 3 * y }, letterR)
plot(Z)
plot(letterR, add=TRUE)

plot(blur(Z, 0.3, bleed=TRUE))
plot(letterR, add=TRUE)

plot(blur(Z, 0.3, bleed=FALSE))
plot(letterR, add=TRUE)
          
plot(blur(Z, 0.3, bleed=FALSE))
plot(letterR, add=TRUE)
          
par(oldpar)

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
applynbd(redwood, R=0.2, showoffK, fullpicture=redwood, rad=0.2, exclude=TRUE)

options(oldoptions)

