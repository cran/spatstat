if(dev.cur() <= 1) get(getOption("device"))()

oldpar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "windows", "Macintosh")))

data(swedishpines)
plot(swedishpines, main="Swedish Pines")

K <- Kest(swedishpines)

conspire(K, cbind(border,theo)~r, subset="r <= 40")
title(main="K function for Swedish Pines")

e <- exactdt(swedishpines)
image(e$xcol, e$yrow, t(e$d), axes=F, xlab="", ylab="",
      main="Distance transform")
points(swedishpines)

image(e$xcol, e$yrow, t(e$d < 4.5), axes=F, xlab="", ylab="")
points(swedishpines)

plot(allstats(swedishpines),
     subset=list("r <= 15", "r <= 15", "r <= 9", "r <= 30"))

fit <- mpl(swedishpines, ~1, Strauss(r=7))
print(fit)

Xsim <- rmh("strauss", c(0.03,0.2,7), swedishpines$window,
            n.start=swedishpines$n, nrep=1e4)

plot(Xsim, main="Simulation from fitted Strauss model")

data(demopat)
plot(demopat, box=FALSE)
plot(alltypes(demopat, "K"), subset="r <= 1500")

fit <- mpl(demopat, ~marks + polynom(x,y,2), Poisson())
plot(fit)

par(oldpar)

