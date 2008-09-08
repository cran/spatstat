if(dev.cur() <= 1) get(getOption("device"))()

oldpar <- par(ask = interactive() && dev.interactive(orNone=TRUE))
oldoptions <- options(warn=-1)

data(amacrine)
plot(amacrine)

data(anemones)
plot(anemones, markscale=0.5)

data(ants)
ants.extra$plot()

data(bei)
plot(bei.extra$elev, main="Beilschmiedia")
plot(bei, add=TRUE, pch=16, cex=0.3)

data(betacells)
plot(betacells)

data(bramblecanes)
plot(bramblecanes, cols=1:3)
plot(split(bramblecanes))

data(cells)
plot(cells)

data(chorley)
chorley.extra$plotit()

data(copper)
plot(copper$Points, main="Copper")
plot(copper$Lines, add=TRUE)

data(demopat)
plot(demopat)

data(finpines)
plot(finpines %mark% finpines.extra$diameter, main="Finnish pines")

data(hamster)
plot(hamster)

data(humberside)
plot(humberside)

data(japanesepines)
plot(japanesepines)

data(lansing)
plot(lansing)
plot(split(lansing))

data(longleaf)
plot(longleaf)

data(murchison)
plot(murchison$greenstone, main="Murchison data", col="lightgreen")
plot(murchison$gold, add=TRUE, pch="+",col="blue")
plot(murchison$faults, add=TRUE, col="red")

data(nbfires)
plot(split(nbfires)$nbfires.00 %mark% factor(nbextras$extras.00$fire.type),
     main=c("New Brunswick fires 2000", "by fire type"),
     cols=c("red", "blue", "green", "cyan"))

data(nztrees)
plot(nztrees)

data(ponderosa)
ponderosa.extra$plotit()

data(redwood)
plot(redwood)

data(redwoodfull)
redwoodfull.extra$plot()

data(residualspaper)
plot(residualspaper$Fig1)
plot(residualspaper$Fig4a)
plot(residualspaper$Fig4b)
plot(residualspaper$Fig4c)


data(simdat)
plot(simdat)

data(spruces)
plot(spruces, maxsize=min(nndist(spruces))/2)

data(swedishpines)
plot(swedishpines)

data(urkiola)
plot(urkiola, cex=0.5)

par(oldpar)
options(oldoptions)
