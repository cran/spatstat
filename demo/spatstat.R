if(dev.cur() <= 1) get(getOption("device"))()

oldpar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "windows", "Macintosh")))
oldoptions <- options(warn=-1)

data(swedishpines)
plot(swedishpines, main="Point pattern")
data(demopat)
plot(demopat, cols=c("green", "blue"), main="Multitype point pattern")
data(longleaf)
plot(longleaf, fg="blue", main="Marked point pattern")

data(lansing)
plot(lansing, "Lansing Woods data")
plot(split(lansing))

data(letterR)
plot(letterR)
lambda <- 10/area.owin(letterR)
points(rpoispp(lambda, win=letterR))
points(rpoispp(9 * lambda, win=letterR))
points(rpoispp(90 * lambda, win=letterR))

X <- swedishpines
subset <- 1:20
plot(X[subset])
subwindow <- owin(poly=list(x=c(0,96,96,40,40,0),y=c(0,0,100,100,50,0)))
plot(X[,subwindow])

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

plot(allstats(swedishpines))

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

par(oldpar)

showoffK <- function(Y, u, d, r, fullpicture,rad) { 
	plot(fullpicture,
             main="Animation using \`applynbd\'\n explaining the K function")
	points(Y[r>0], cex=2)
	points(u[1],u[2],pch="+",cex=3)
	theta <- seq(0,2*pi,length=100)
	polygon(u[1]+ rad * cos(theta),u[2]+rad*sin(theta))
	text(u[1]+rad/3,u[2]+rad/2,Y$n-1,cex=3)
	Sys.sleep(if(runif(1) < 0.05) 1.2 else 0.25)
	return(Y$n - 1)
}
data(redwood)
applynbd(redwood, R=0.2, showoffK, fullpicture=redwood, rad=0.2)

options(oldoptions)
