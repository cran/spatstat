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

X <- swedishpines
subset <- 1:20
plot(X[subset])
subwindow <- owin(poly=list(x=c(0,96,96,40,40,0),y=c(0,0,100,100,50,0)))
plot(X[,subwindow])

K <- Kest(swedishpines)
conspire(K, cbind(border,theo)~r, subset="r <= 40")
title(main="K function for Swedish Pines")

e <- exactdt(swedishpines)
image(e$xcol, e$yrow, t(e$d), axes=FALSE, xlab="", ylab="",
      main="Distance transform")
points(swedishpines)

image(e$xcol, e$yrow, t(e$d < 4.5), axes=FALSE, xlab="", ylab="")
points(swedishpines)

plot(allstats(swedishpines),
     subset=list("r <= 15", "r <= 15", "r <= 9", "r <= 30"))

fit <- mpl(swedishpines, ~1, Strauss(r=7))
print(fit)

Xsim <- rmh("strauss", c(0.03,0.2,7), swedishpines$window,
            n.start=swedishpines$n, nrep=1e4)

plot(Xsim, main="Simulation from fitted Strauss model")

data(demopat)
plot(demopat, cols=c("red", "blue"))
plot(alltypes(demopat, "K"), subset="r <= 1500")

fit <- mpl(demopat, ~marks + polynom(x,y,2), Poisson())
plot(fit)

plot(rpoispp(100))
plot(rpoispp(function(x,y){1000 * exp(-3*x)}, 1000))
plot(rMaternII(200, 0.05))
plot(rSSI(0.05, 200))
plot(rThomas(10, 0.2, 5))
plot(rMatClust(10, 0.05, 4))
Xg <- rmh("geyer", par=c(1.25, 1.6, 0.2, 4.5),w=c(0,10,0,10),
          n.start=200, nrep=1e4)
plot(Xg, main="rmh(\"geyer\", ...)")

par(oldpar)
options(oldoptions)
