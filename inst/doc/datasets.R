### R code from vignette source 'datasets.Rnw'

###################################################
### code chunk number 1: datasets.Rnw:5-6
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: datasets.Rnw:27-34
###################################################
library(spatstat)
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
spatstat.options(transparent=FALSE)
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: datasets.Rnw:223-241
###################################################
opa <- par()
## How to set all margins to zero and eliminate all outer spaces
zeromargins <- function() {
  par(
      mar=rep(0,4),
      omd=c(0,1,0,1),
      xaxs="i",
      yaxs="i"
  )
  invisible(NULL)
}
## Set 'mar'
setmargins <- function(...) {
  x <- c(...)
  x <- rep(x, 4)[1:4]
  par(mar=x)
  invisible(NULL)
}


###################################################
### code chunk number 4: datasets.Rnw:250-251 (eval = FALSE)
###################################################
## plot(amacrine)


###################################################
### code chunk number 5: datasets.Rnw:253-255
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,1,2,0)
plot(amacrine)


###################################################
### code chunk number 6: datasets.Rnw:264-265 (eval = FALSE)
###################################################
## plot(anemones, markscale=1)


###################################################
### code chunk number 7: datasets.Rnw:267-269
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,0,2,0)
plot(anemones, markscale=1)


###################################################
### code chunk number 8: datasets.Rnw:282-283 (eval = FALSE)
###################################################
## ants.extra$plotit()


###################################################
### code chunk number 9: datasets.Rnw:285-287
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,0,1,0)
ants.extra$plotit()


###################################################
### code chunk number 10: datasets.Rnw:295-296
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(austates)


###################################################
### code chunk number 11: datasets.Rnw:306-308 (eval = FALSE)
###################################################
## plot(bdspots, equal.scales=TRUE, pch="+", 
##      panel.args=function(i)list(cex=c(0.15, 0.2, 0.7)[i]))


###################################################
### code chunk number 12: datasets.Rnw:310-314
###################################################
getOption("SweaveHooks")[["fig"]]()
zeromargins()
plot(bdspots, equal.scales=TRUE, pch="+", main="",
     mar.panel=0, hsep=1,
     panel.args=function(i)list(cex=c(0.15, 0.2, 0.7)[i]))


###################################################
### code chunk number 13: datasets.Rnw:324-326 (eval = FALSE)
###################################################
## plot(bei.extra$elev, main="Beilschmiedia")
## plot(bei, add=TRUE, pch=16, cex=0.3)


###################################################
### code chunk number 14: datasets.Rnw:328-331
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,0,2,0)
plot(bei.extra$elev, main="Beilschmiedia")
plot(bei, add=TRUE, pch=16, cex=0.3)


###################################################
### code chunk number 15: datasets.Rnw:337-343 (eval = FALSE)
###################################################
## M <- persp(bei.extra$elev, 
##            theta=-45, phi=18, expand=7,
##            border=NA, apron=TRUE, shade=0.3, 
##            box=FALSE, visible=TRUE,
##            main="")
## perspPoints(bei, Z=bei.extra$elev, M=M, pch=16, cex=0.3)


###################################################
### code chunk number 16: datasets.Rnw:352-353
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(betacells)


###################################################
### code chunk number 17: datasets.Rnw:358-359
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(bramblecanes, cols=1:3)


###################################################
### code chunk number 18: datasets.Rnw:364-365 (eval = FALSE)
###################################################
## plot(split(bramblecanes))


###################################################
### code chunk number 19: datasets.Rnw:375-376
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(bronzefilter,markscale=2)


###################################################
### code chunk number 20: datasets.Rnw:384-385
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(btb, which.marks="spoligotype", cols=2:5, chars=1:4)


###################################################
### code chunk number 21: datasets.Rnw:394-395
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(cells)


###################################################
### code chunk number 22: datasets.Rnw:403-404
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(cetaceans.extra$patterns, main="Cetaceans data", cols=1:5, hsep=1)


###################################################
### code chunk number 23: datasets.Rnw:413-416
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(chicago, main="Chicago Crimes", col="grey",
     cols=c("red", "blue", "black", "blue", "red", "blue", "blue"),
     chars=c(16,2,22,17,24,15,6), leg.side="left", show.window=FALSE)


###################################################
### code chunk number 24: datasets.Rnw:426-427
###################################################
getOption("SweaveHooks")[["fig"]]()
chorley.extra$plotit()


###################################################
### code chunk number 25: datasets.Rnw:443-445
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(clmfires, which.marks="cause", cols=2:5, cex=0.25,
     main="Castilla-La Mancha forest fires")


###################################################
### code chunk number 26: datasets.Rnw:455-456
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(clmfires.extra$clmcov100$elevation, main="Elevation")


###################################################
### code chunk number 27: datasets.Rnw:467-468
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(concrete,chars="+",cols="blue",col="yellow")


###################################################
### code chunk number 28: datasets.Rnw:479-481
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(copper$Points, main="Copper")
plot(copper$Lines, add=TRUE)


###################################################
### code chunk number 29: datasets.Rnw:488-490
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }),
      parargs=list(mar=rep(1,4)))


###################################################
### code chunk number 30: datasets.Rnw:497-498
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demopat)


###################################################
### code chunk number 31: datasets.Rnw:512-513
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(dendrite, leg.side="bottom", main="", cex=0.75, cols=2:4)


###################################################
### code chunk number 32: datasets.Rnw:521-522
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(finpines, main="Finnish pines")


###################################################
### code chunk number 33: datasets.Rnw:535-539
###################################################
getOption("SweaveHooks")[["fig"]]()
wildM1 <- with(flu, virustype == "wt" & stain == "M2-M1")
plot(flu[wildM1, 1, drop=TRUE],
     main=c("flu data", "wild type virus, M2-M1 stain"),
     chars=c(16,3), cex=0.4, cols=2:3)


###################################################
### code chunk number 34: datasets.Rnw:547-548
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(gordon, main="People in Gordon Square", pch=16)


###################################################
### code chunk number 35: datasets.Rnw:563-564
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(gorillas, which.marks=1, chars=c(1,3), cols=2:3, main="Gorilla nest sites")


###################################################
### code chunk number 36: datasets.Rnw:568-569 (eval = FALSE)
###################################################
## system.file("rawdata/gorillas/vegetation.asc", package="spatstat")


###################################################
### code chunk number 37: datasets.Rnw:578-579
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(hamster, cols=c(2,4))


###################################################
### code chunk number 38: datasets.Rnw:589-590
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(heather$coarse)


###################################################
### code chunk number 39: datasets.Rnw:594-595 (eval = FALSE)
###################################################
## plot(heather)


###################################################
### code chunk number 40: datasets.Rnw:605-606
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(humberside)


###################################################
### code chunk number 41: datasets.Rnw:618-619
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(hyytiala, cols=2:5)


###################################################
### code chunk number 42: datasets.Rnw:628-629
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(japanesepines)


###################################################
### code chunk number 43: datasets.Rnw:638-639
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(lansing)


###################################################
### code chunk number 44: datasets.Rnw:645-646 (eval = FALSE)
###################################################
## plot(split(lansing))


###################################################
### code chunk number 45: datasets.Rnw:653-654
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(longleaf)


###################################################
### code chunk number 46: datasets.Rnw:663-668
###################################################
getOption("SweaveHooks")[["fig"]]()
pa <- function(i) {
   if(i == 1) list(cols=c("red", "green")) else 
   list(do.col=TRUE, col=grey(seq(1,0,length=32)))
}
plot(meningitis, panel.args=pa)


###################################################
### code chunk number 47: datasets.Rnw:677-679
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mucosa, chars=c(1,3), cols=c("red", "green"))
plot(mucosa.subwin, add=TRUE, lty=3)


###################################################
### code chunk number 48: datasets.Rnw:695-698
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(murchison$greenstone, main="Murchison data", col="lightgreen")
plot(murchison$gold, add=TRUE, pch=3, col="blue")
plot(murchison$faults, add=TRUE, col="red")


###################################################
### code chunk number 49: datasets.Rnw:704-705
###################################################
reedy <- owin(c(580, 650) * 1000, c(6986, 7026) * 1000)


###################################################
### code chunk number 50: datasets.Rnw:710-713
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(murchison$greenstone[reedy], main="Murchison data", col="lightgreen")
plot(murchison$gold[reedy], add=TRUE, pch=3, col="blue")
plot(murchison$faults[reedy], add=TRUE, col="red")


###################################################
### code chunk number 51: datasets.Rnw:721-722
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(nbfires, use.marks=FALSE, pch=".")


###################################################
### code chunk number 52: datasets.Rnw:728-729 (eval = FALSE)
###################################################
## plot(split(nbfires), use.marks=FALSE, chars=".")


###################################################
### code chunk number 53: datasets.Rnw:732-737
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(0,0,2,0))
plot(split(nbfires)$"2000", which.marks="fire.type",
     main=c("New Brunswick fires 2000", "by fire type"),
     cols=c("blue", "green", "red", "cyan"),
     leg.side="left")


###################################################
### code chunk number 54: datasets.Rnw:745-747
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(nztrees)
plot(trim.rectangle(as.owin(nztrees), c(0,5), 0), add=TRUE, lty=3)


###################################################
### code chunk number 55: datasets.Rnw:760-761
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(osteo[1:4,], main.panel="", pch=21, bg='white')


###################################################
### code chunk number 56: datasets.Rnw:767-768 (eval = FALSE)
###################################################
## system.file("rawdata/osteo/osteo36.txt", package="spatstat")


###################################################
### code chunk number 57: datasets.Rnw:777-778
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(paracou, cols=2:3, chars=c(16,3))


###################################################
### code chunk number 58: datasets.Rnw:786-787
###################################################
getOption("SweaveHooks")[["fig"]]()
ponderosa.extra$plotit()


###################################################
### code chunk number 59: datasets.Rnw:799-800
###################################################
pyr <- pyramidal[c(FALSE,TRUE), ]


###################################################
### code chunk number 60: datasets.Rnw:803-805
###################################################
getOption("SweaveHooks")[["fig"]]()
pyr$grp <- abbreviate(pyr$group, minlength=7)
plot(pyr, quote(plot(Neurons, pch=16, main=grp)), main="Pyramidal Neurons")


###################################################
### code chunk number 61: datasets.Rnw:825-827
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(redwood)
plot(redwood3, add=TRUE, pch=20)


###################################################
### code chunk number 62: datasets.Rnw:830-831
###################################################
getOption("SweaveHooks")[["fig"]]()
redwoodfull.extra$plotit()


###################################################
### code chunk number 63: datasets.Rnw:845-847
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(as.solist(residualspaper[c("Fig1", "Fig4a", "Fig4b", "Fig4c")]), 
     main="")


###################################################
### code chunk number 64: datasets.Rnw:855-856
###################################################
getOption("SweaveHooks")[["fig"]]()
shapley.extra$plotit(main="Shapley")


###################################################
### code chunk number 65: datasets.Rnw:863-865
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(shelling, pch=3)
plot(onearrow(830, 400, 830, 530, "N"), add=TRUE)


###################################################
### code chunk number 66: datasets.Rnw:872-873
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simdat)


###################################################
### code chunk number 67: datasets.Rnw:881-882
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(spiders, pch=16, show.window=FALSE)


###################################################
### code chunk number 68: datasets.Rnw:889-892
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(sporophores, chars=c(16,1,2), cex=0.6)
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)


###################################################
### code chunk number 69: datasets.Rnw:904-905
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(stonetools, which.marks=2, cols=c(2,3), chars=c(1,3), cex=0.5)


###################################################
### code chunk number 70: datasets.Rnw:914-915
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(spruces, maxsize=min(nndist(spruces)))


###################################################
### code chunk number 71: datasets.Rnw:924-925
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(swedishpines)


###################################################
### code chunk number 72: datasets.Rnw:934-935
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(urkiola, cex=0.5, cols=2:3)


###################################################
### code chunk number 73: datasets.Rnw:942-944
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(0,0,2,0))
plot(waka, markscale=0.04, main=c("Waka national park", "tree diameters"))


###################################################
### code chunk number 74: datasets.Rnw:951-955
###################################################
getOption("SweaveHooks")[["fig"]]()
v <- rotate(vesicles, pi/2)
ve <- lapply(vesicles.extra, rotate, pi/2)
plot(v, main="Vesicles")
plot(ve$activezone, add=TRUE, lwd=3)


###################################################
### code chunk number 75: datasets.Rnw:980-981 (eval = FALSE)
###################################################
## system.file("rawdata/vesicles/mitochondria.txt", package="spatstat")


###################################################
### code chunk number 76: datasets.Rnw:989-990
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders)


