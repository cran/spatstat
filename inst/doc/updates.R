### R code from vignette source 'updates.Rnw'

###################################################
### code chunk number 1: updates.Rnw:20-24
###################################################
library(spatstat)
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: updates.Rnw:38-117
###################################################
readSizeTable <- function(fname) {
  if(is.null(fname) || !file.exists(fname)) return(NULL)
  a <- read.table(fname, header=TRUE)
  a$date <- as.Date(a$date)
  return(a)
}
getSizeTable <- function(packagename="spatstat", tablename="packagesizes.txt") {
  fname <- system.file("doc", tablename, package=packagename)
  readSizeTable(fname)
}
counts <- c("nhelpfiles", "nobjects", "ndatasets", "Rlines", "srclines")
mergeSizeTables <- function(a, b, breakupdate) {
  #' a is the running total for spatstat; b is a sub-package.
  #' breakupdate is the date when the code in b was removed from spatstat
  #' so that the size of 'b' must be added to 'a' for all dates >= breakupdate
  if(is.null(b)) return(a)
  adates <- a$date
  bdates <- b$date
  alldates <- sort(unique(c(adates,bdates)))
  if(missing(breakupdate)) breakupdate <- min(bdates)
  #' functions to determine, for any given date,
  #' the relevant (latest) row of the table
  aok <- rev(!duplicated(rev(adates)))
  arowfun <- approxfun(adates[aok], seq_along(adates)[aok], 
                       method="constant", f=0, rule=2, yleft=0)
  bok <- rev(!duplicated(rev(bdates)))
  browfun <- approxfun(bdates[bok], seq_along(bdates)[bok], 
                       method="constant", f=0, rule=2, yleft=0)
  result <- NULL
  for(k in seq_along(alldates)) {
    thedate <- alldates[k]
    i <- arowfun(thedate)
    j <- browfun(thedate)
    #' i > 0 because spatstat's founding date is earlier than any sub-package
    nextrow <- a[i, ]
    if(j > 0 && thedate >= breakupdate) {
      #' add contribution from 'b'
      nextrow[, counts] <- nextrow[, counts] + b[j, counts]
    }
    result <- rbind(result, nextrow)
  }
  return(result)
}
z <- getSizeTable()
## handle the case where spatstat is the umbrella package
isumb <- existsSpatstatVariable("Spatstat.Is.Umbrella") && 
         isTRUE(getSpatstatVariable("Spatstat.Is.Umbrella"))
if(isumb) {
  Bday <- "2020-12-14"
  zgeom <- getSizeTable("spatstat.geom")
  z <- mergeSizeTables(z, zgeom, Bday)
  zcore <- getSizeTable("spatstat.core")
  z <- mergeSizeTables(z, zcore, Bday)
  zlin <- getSizeTable("spatstat.linnet")
  z <- mergeSizeTables(z, zlin, Bday)
} 
## sub-packages- access via the sub-package
zutils <- getSizeTable("spatstat.utils")
zdata <- getSizeTable("spatstat.data")
zsparse <- getSizeTable("spatstat.sparse")
z <- mergeSizeTables(z, zutils, "2017-03-22")
z <- mergeSizeTables(z, zdata,  "2017-09-23")
z <- mergeSizeTables(z, zsparse, "2020-11-04")
## extension packages - use copy of package size file stored in spatstat
zlocal <- getSizeTable("spatstat", "spatstatlocalsize.txt")
zgui <- getSizeTable("spatstat", "spatstatguisize.txt")
zKnet <- getSizeTable("spatstat", "spatstatKnetsize.txt")
z <- mergeSizeTables(z, zlocal)
z <- mergeSizeTables(z, zgui)
z <- mergeSizeTables(z, zKnet)
#
currentcount <- z[nrow(z), counts]
bookcount    <- z[z$version == "1.42-0", counts]
changes <- currentcount - bookcount
newobj <- changes[["nobjects"]]
newdat <- changes[["ndatasets"]] + 1  # counting rule doesn't detect redwood3
newcode  <- changes[["Rlines"]] + changes[["srclines"]]
bookcode <- bookcount[["Rlines"]] + bookcount[["srclines"]]
growth <- signif((100 * newcode)/bookcode, digits=2)


###################################################
### code chunk number 3: updates.Rnw:134-139
###################################################
options(SweaveHooks=list(fig=function() par(mar=0.2+c(2,4,2,0))))
Plot <- function(fmla, ..., dat=z) {
  yvals <- eval(as.expression(fmla[[2]]), envir=dat)
  plot(fmla, ..., data=dat, type="l", xlab="", lwd=2, ylim=c(0, max(yvals)))
}


###################################################
### code chunk number 4: updates.Rnw:145-150
###################################################
getOption("SweaveHooks")[["fig"]]()
Plot((Rlines + srclines)/1000 ~ date, ylab="Lines of code (x 1000)", 
     main="Spatstat growth")
lines(srclines/1000 ~ date, data=z)
text(as.Date("2015-01-01"), 9.5, "C code")
text(as.Date("2015-01-01"), 60, "R code")


