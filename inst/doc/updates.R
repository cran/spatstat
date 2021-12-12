### R code from vignette source 'updates.Rnw'

###################################################
### code chunk number 1: updates.Rnw:22-28
###################################################
library(spatstat)
x <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
              fields = c("Version", "Date"))
sversion <- as.character(x[,"Version"])
sdate    <- as.character(x[,"Date"])
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: updates.Rnw:38-125
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
RemoveDevel <- function(sizetable) {
  ## remove entries with fractional version numbers
  if(is.null(sizetable)) return(NULL)
  ver <- sizetable$version
  isdevel <- sapply(ver, function(x) { length(unlist(package_version(x))) > 3 })
  st <- if(all(isdevel)) NULL else sizetable[!isdevel, , drop=FALSE]
  return(st)
}
counts <- c("nhelpfiles", "nobjects", "ndatasets", "Rlines", "srclines")
mergeSizeTables <- function(a, b, breakupdate, allow.devel=FALSE) {
  #' a is the running total for spatstat; b is a sub-package.
  #' breakupdate is the date when the code in b was removed from spatstat
  #' so that the size of 'b' must be added to 'a' for all dates >= breakupdate
  if(!allow.devel) b <- RemoveDevel(b)
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
## start collating size tables for the spatstat family
z <- getSizeTable()
## sub-packages- access via the sub-package
zgeom <- getSizeTable("spatstat.geom")
zcore <- getSizeTable("spatstat.core")
zlin <- getSizeTable("spatstat.linnet")
zutils <- getSizeTable("spatstat.utils")
zdata <- getSizeTable("spatstat.data")
zsparse <- getSizeTable("spatstat.sparse")
## merge history of geom, core, linnet after the big split
Bday <- "2020-12-14"
z <- mergeSizeTables(z, zgeom, Bday)
z <- mergeSizeTables(z, zcore, Bday)
z <- mergeSizeTables(z, zlin, Bday)
## merge history of other packages after their split dates
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
currentcode <- currentcount[["Rlines"]] + currentcount[["srclines"]]
growth <- signif((100 * newcode)/bookcode, digits=2)


###################################################
### code chunk number 3: updates.Rnw:137-142
###################################################
options(SweaveHooks=list(fig=function() par(mar=0.2+c(2,4,2,0))))
Plot <- function(fmla, ..., dat=z) {
  yvals <- eval(as.expression(fmla[[2]]), envir=dat)
  plot(fmla, ..., data=dat, type="l", xlab="", lwd=2, ylim=c(0, max(yvals)))
}


###################################################
### code chunk number 4: updates.Rnw:148-153
###################################################
getOption("SweaveHooks")[["fig"]]()
Plot((Rlines + srclines)/1000 ~ date, ylab="Lines of code (x 1000)", 
     main="Spatstat growth")
lines(srclines/1000 ~ date, data=z)
text(as.Date("2015-01-01"), 9.5, "C code")
text(as.Date("2015-01-01"), 60, "R code")


###################################################
### code chunk number 5: updates.Rnw:170-188
###################################################
## Tabulate latest version numbers of packages
vtable <- data.frame(package="spatstat", version=sversion, date=as.Date(sdate))
AppendVersion <- function(pkg, sizetable, v, allow.devel=FALSE) {
  if(!allow.devel) sizetable <- RemoveDevel(sizetable)
  if(is.null(sizetable)) return(v)
  lastrow <- sizetable[nrow(sizetable), , drop=FALSE]
  if(is.null(lastrow)) return(v)
  rbind(v, data.frame(package=pkg, version=lastrow[,"version"], date=as.Date(lastrow[,"date"])))
}
vtable <- AppendVersion("spatstat.geom", zgeom, vtable)
vtable <- AppendVersion("spatstat.core", zcore, vtable)
vtable <- AppendVersion("spatstat.linnet", zlin, vtable)
vtable <- AppendVersion("spatstat.sparse", zsparse, vtable)
vtable <- AppendVersion("spatstat.data", zdata, vtable)
vtable <- AppendVersion("spatstat.utils", zutils, vtable)
vtable <- AppendVersion("spatstat.local", zlocal, vtable)
vtable <- AppendVersion("spatstat.Knet", zKnet, vtable)
vtable <- AppendVersion("spatstat.gui", zgui, vtable)


###################################################
### code chunk number 6: updates.Rnw:194-195
###################################################
print(vtable[,c(3,1,2)], row.names=FALSE)


