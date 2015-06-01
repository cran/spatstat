### R code from vignette source 'datasets.Rnw'

###################################################
### code chunk number 1: datasets.Rnw:5-6
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: datasets.Rnw:25-31
###################################################
library(spatstat)
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: datasets.Rnw:199-218
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
  par(opa)
  x <- c(...)
  x <- rep(x, 4)[1:4]
  par(mar=x)
  invisible(NULL)
}


###################################################
### code chunk number 4: datasets.Rnw:227-228 (eval = FALSE)
###################################################
## plot(amacrine)


###################################################
### code chunk number 5: datasets.Rnw:230-232
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,1,0,0)
plot(amacrine)


###################################################
### code chunk number 6: datasets.Rnw:241-242 (eval = FALSE)
###################################################
## plot(anemones, markscale=1)


###################################################
### code chunk number 7: datasets.Rnw:244-246
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0)
plot(anemones, markscale=1)


###################################################
### code chunk number 8: datasets.Rnw:259-260 (eval = FALSE)
###################################################
## ants.extra$plotit()


###################################################
### code chunk number 9: datasets.Rnw:262-264
###################################################
getOption("SweaveHooks")[["fig"]]()
zeromargins()
ants.extra$plotit()


###################################################
### code chunk number 10: datasets.Rnw:274-276 (eval = FALSE)
###################################################
## plot(bdspots, equal.scales=TRUE, pch="+", 
##      panel.args=function(i)list(cex=c(0.15, 0.2, 0.7)[i]))


###################################################
### code chunk number 11: datasets.Rnw:278-282
###################################################
getOption("SweaveHooks")[["fig"]]()
zeromargins()
plot(bdspots, equal.scales=TRUE, pch="+", main="",
     mar.panel=0, hsep=1,
     panel.args=function(i)list(cex=c(0.15, 0.2, 0.7)[i]))


###################################################
### code chunk number 12: datasets.Rnw:292-294 (eval = FALSE)
###################################################
## plot(bei.extra$elev, main="Beilschmiedia")
## plot(bei, add=TRUE, pch=16, cex=0.3)


###################################################
### code chunk number 13: datasets.Rnw:296-299
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0)
plot(bei.extra$elev, main="Beilschmiedia")
plot(bei, add=TRUE, pch=16, cex=0.3)


###################################################
### code chunk number 14: datasets.Rnw:302-308
###################################################
getOption("SweaveHooks")[["fig"]]()
M <- persp(bei.extra$elev, 
           theta=-45, phi=18, expand=7,
           border=NA, apron=TRUE, shade=0.3, 
           box=FALSE, visible=TRUE,
           main="")
perspPoints(bei, Z=bei.extra$elev, M=M, pch=16, cex=0.3)


###################################################
### code chunk number 15: datasets.Rnw:317-318
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(betacells)


###################################################
### code chunk number 16: datasets.Rnw:323-324
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(bramblecanes, cols=1:3)


###################################################
### code chunk number 17: datasets.Rnw:327-328
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(bramblecanes))


###################################################
### code chunk number 18: datasets.Rnw:338-339
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(bronzefilter,markscale=2)


###################################################
### code chunk number 19: datasets.Rnw:348-349
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(cells)


###################################################
### code chunk number 20: datasets.Rnw:358-361
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(chicago, main="Chicago Street Crimes", col="grey",
     cols=c("red", "blue", "black", "blue", "red", "blue", "blue"),
     chars=c(16,2,22,17,24,15,6), leg.side="left", show.window=FALSE)


###################################################
### code chunk number 21: datasets.Rnw:371-372
###################################################
getOption("SweaveHooks")[["fig"]]()
chorley.extra$plotit()


###################################################
### code chunk number 22: datasets.Rnw:388-390
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(clmfires, which.marks="cause", cols=2:5, cex=0.25,
     main="Castilla-La Mancha forest fires")


###################################################
### code chunk number 23: datasets.Rnw:400-401
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(clmfires.extra$clmcov200, main="Covariates for forest fires")


###################################################
### code chunk number 24: datasets.Rnw:412-414
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(copper$Points, main="Copper")
plot(copper$Lines, add=TRUE)


###################################################
### code chunk number 25: datasets.Rnw:421-423
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }),
      parargs=list(mar=rep(1,4)))


###################################################
### code chunk number 26: datasets.Rnw:430-431
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demopat)


###################################################
### code chunk number 27: datasets.Rnw:445-446
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(dendrite, leg.side="bottom", main="", cex=0.75, cols=2:4)


###################################################
### code chunk number 28: datasets.Rnw:454-455
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(finpines, main="Finnish pines")


###################################################
### code chunk number 29: datasets.Rnw:468-472
###################################################
getOption("SweaveHooks")[["fig"]]()
wildM1 <- with(flu, virustype == "wt" & stain == "M2-M1")
plot(flu[wildM1, 1, drop=TRUE],
     main=c("flu data", "wild type virus, M2-M1 stain"),
     chars=c(16,3), cex=0.4, cols=2:3)


###################################################
### code chunk number 30: datasets.Rnw:480-481
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(gordon, main="People in Gordon Square", pch=16)


###################################################
### code chunk number 31: datasets.Rnw:496-497
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(gorillas, which.marks=1, chars=c(1,3), cols=2:3, main="Gorilla nest sites")


###################################################
### code chunk number 32: datasets.Rnw:501-502 (eval = FALSE)
###################################################
## system.file("rawdata/gorillas/vegetation.asc", package="spatstat")


###################################################
### code chunk number 33: datasets.Rnw:511-512
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(hamster, cols=c(2,4))


###################################################
### code chunk number 34: datasets.Rnw:522-523
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(heather)


###################################################
### code chunk number 35: datasets.Rnw:533-534
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(humberside)


###################################################
### code chunk number 36: datasets.Rnw:546-547
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(hyytiala, cols=2:5)


###################################################
### code chunk number 37: datasets.Rnw:556-557
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(japanesepines)


###################################################
### code chunk number 38: datasets.Rnw:566-567
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(lansing)


###################################################
### code chunk number 39: datasets.Rnw:570-571
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(lansing))


###################################################
### code chunk number 40: datasets.Rnw:578-579
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(longleaf)


###################################################
### code chunk number 41: datasets.Rnw:588-590
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mucosa, chars=c(1,3), cols=c("red", "green"))
plot(mucosa.subwin, add=TRUE, lty=3)


###################################################
### code chunk number 42: datasets.Rnw:604-607
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(murchison$greenstone, main="Murchison data", col="lightgreen")
plot(murchison$gold, add=TRUE, pch=3, col="blue")
plot(murchison$faults, add=TRUE, col="red")


###################################################
### code chunk number 43: datasets.Rnw:615-616
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(nbfires, use.marks=FALSE, pch=".")


###################################################
### code chunk number 44: datasets.Rnw:619-620
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(nbfires), use.marks=FALSE, chars=".")


###################################################
### code chunk number 45: datasets.Rnw:623-627
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(nbfires)$"2000", which.marks="fire.type",
     main=c("New Brunswick fires 2000", "by fire type"),
     cols=c("blue", "green", "red", "cyan"),
     leg.side="left")


###################################################
### code chunk number 46: datasets.Rnw:635-637
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(nztrees)
plot(trim.rectangle(as.owin(nztrees), c(0,5), 0), add=TRUE, lty=3)


###################################################
### code chunk number 47: datasets.Rnw:650-651
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(osteo[1:10,], main.panel="", pch=21, bg='white')


###################################################
### code chunk number 48: datasets.Rnw:657-658 (eval = FALSE)
###################################################
## system.file("rawdata/osteo/osteo36.txt", package="spatstat")


###################################################
### code chunk number 49: datasets.Rnw:667-668
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(paracou, cols=2:3, chars=c(16,3))


###################################################
### code chunk number 50: datasets.Rnw:676-677
###################################################
getOption("SweaveHooks")[["fig"]]()
ponderosa.extra$plotit()


###################################################
### code chunk number 51: datasets.Rnw:688-691
###################################################
getOption("SweaveHooks")[["fig"]]()
pyr <- pyramidal
pyr$grp <- abbreviate(pyramidal$group, minlength=7)
plot(pyr, quote(plot(Neurons, pch=16, main=grp)), main="Pyramidal Neurons")


###################################################
### code chunk number 52: datasets.Rnw:708-709
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(redwood)


###################################################
### code chunk number 53: datasets.Rnw:712-713
###################################################
getOption("SweaveHooks")[["fig"]]()
redwoodfull.extra$plotit()


###################################################
### code chunk number 54: datasets.Rnw:727-729
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(as.listof(residualspaper[c("Fig1", "Fig4a", "Fig4b", "Fig4c")]), 
     main="")


###################################################
### code chunk number 55: datasets.Rnw:737-738
###################################################
getOption("SweaveHooks")[["fig"]]()
shapley.extra$plotit(main="Shapley")


###################################################
### code chunk number 56: datasets.Rnw:745-746
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simdat)


###################################################
### code chunk number 57: datasets.Rnw:754-755
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(spiders, pch=16, show.window=FALSE)


###################################################
### code chunk number 58: datasets.Rnw:762-765
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(sporophores, chars=c(16,1,2), cex=0.6)
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)


###################################################
### code chunk number 59: datasets.Rnw:774-775
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(spruces, maxsize=min(nndist(spruces)))


###################################################
### code chunk number 60: datasets.Rnw:784-785
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(swedishpines)


###################################################
### code chunk number 61: datasets.Rnw:794-795
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(urkiola, cex=0.5, cols=2:3)


###################################################
### code chunk number 62: datasets.Rnw:802-803
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waka, markscale=0.04, main=c("Waka national park", "tree diameters"))


###################################################
### code chunk number 63: datasets.Rnw:810-814
###################################################
getOption("SweaveHooks")[["fig"]]()
v <- rotate(vesicles, pi/2)
ve <- lapply(vesicles.extra, rotate, pi/2)
plot(v, main="Vesicles")
plot(ve$activezone, add=TRUE, lwd=3)


###################################################
### code chunk number 64: datasets.Rnw:839-840 (eval = FALSE)
###################################################
## system.file("rawdata/vesicles/mitochondria.txt", package="spatstat")


###################################################
### code chunk number 65: datasets.Rnw:848-849
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders)


