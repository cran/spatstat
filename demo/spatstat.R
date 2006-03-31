if(dev.cur() <= 1) get(getOption("device"))()

oldpar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "windows", "Macintosh")))
oldoptions <- options(warn=-1)

fanfare <- function(stuff) {
  plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
  text(0.5,0.5, stuff, cex=3)
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

Z <- as.im(function(x,y){ sqrt((x - 1)^2 + (y-1)^2)}, square(2))
plot(Z, main="Pixel image")

fanfare("II. Basic operations")
X <- swedishpines
subset <- 1:20
plot(X[subset], main="subset operation: X[subset]")
subwindow <- owin(poly=list(x=c(0,96,96,40,40,0),y=c(0,0,100,100,50,0)))
plot(X[,subwindow], main="subset operation: X[, subwindow]")

data(lansing)
plot(lansing, "Lansing Woods data")
plot(split(lansing))

L <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
w <- owin(c(0.1,0.7), c(0.2, 0.8))
S <- L[, w]
plot(L, main="Subset of a line segment pattern")
plot(S, add=TRUE, col="red")

fanfare("III. Exploratory data analysis")

data(cells)
Z <- density.ppp(cells, 0.07)
plot(Z, main="Kernel smoothed intensity of point pattern")
plot(cells, add=TRUE)

D <- density(L, sigma=0.05)
plot(D, main="Kernel smoothed intensity of line segment pattern")
plot(L, add=TRUE)

plot(swedishpines, main="Swedish Pines data")
K <- Kest(swedishpines)
plot(K)
title(main="K function for Swedish Pines")

en <- envelope(swedishpines, fun=Kest, nsim=10, correction="translate")
plot(en, main="Envelopes of K function based on CSR")

pc <- pcf(swedishpines)
plot(pc)
title(main="Pair correlation function")

plot(swedishpines %mark% (nndist(swedishpines)/2), markscale=1, main="Stienen diagram")

plot(swedishpines$window, main="Distance map")
dis <- distmap(swedishpines)
plot(dis, add=TRUE)
points(swedishpines)

plot(swedishpines$window, main="Thresholded distance")
dis$v <- (dis$v < 4.5)
plot(dis, add=TRUE)
points(swedishpines)

a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
contour(distmap(a), main="Distance map")
plot(a, add=TRUE,col="red")

plot(allstats(swedishpines))

fanfare("IV. Model-fitting")
fit <- ppm(swedishpines, ~1, Strauss(r=7))
print(fit)

Xsim <- rmh(model=fit,
            start=list(n.start=80),
            control=list(nrep=100))
plot(Xsim, main="Simulation from fitted Strauss model")

data(demopat)
plot(demopat, cols=c("red", "blue"))
plot(alltypes(demopat, "K"))

fit <- ppm(demopat, ~marks + polynom(x,y,2), Poisson())
plot(fit)

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
Xg <- rmh(list(cif="geyer", par=c(beta=1.25, gamma=1.6, r=0.2, sat=4.5),
               w=c(0,10,0,10)),
          control=list(nrep=1e4), start=list(n.start=200))
plot(Xg, main=paste("Geyer saturation process\n",
                    "rmh() with cif=\"geyer\""))

plot(rpoisline(10))

fanfare("VI. Programming tools")

par(oldpar)

showoffK <- function(Y, current, ..., fullpicture,rad) { 
	plot(fullpicture,
             main="Animation using \`applynbd\'\n explaining the K function")
	points(Y, cex=2)
        u <- current
	points(u[1],u[2],pch="+",cex=3)
	theta <- seq(0,2*pi,length=100)
	polygon(u[1]+ rad * cos(theta),u[2]+rad*sin(theta))
	text(u[1]+rad/3,u[2]+rad/2,Y$n,cex=3)
        if(runif(1) < 0.2) Sys.sleep(runif(1, max=0.4))
	return(Y$n)
}
data(redwood)
applynbd(redwood, R=0.2, showoffK, fullpicture=redwood, rad=0.2, exclude=TRUE)

options(oldoptions)

