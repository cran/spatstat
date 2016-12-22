### R code from vignette source 'updates.Rnw'

###################################################
### code chunk number 1: updates.Rnw:20-24
###################################################
library(spatstat)
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: updates.Rnw:41-47
###################################################
fname <- system.file("doc", "packagesizes.txt", package="spatstat")
z <- read.table(fname, header=TRUE)
z$date <- as.Date(z$date)
changes <- z[nrow(z), ] - z[z$version == "1.42-0", ]
newobj <- changes[["nobjects"]]
newdat <- changes[["ndatasets"]] + 1  # counting rule doesn't detect redwood3


###################################################
### code chunk number 3: updates.Rnw:57-62
###################################################
options(SweaveHooks=list(fig=function() par(mar=0.2+c(2,4,2,0))))
Plot <- function(fmla, ..., dat=z) {
  yvals <- eval(as.expression(fmla[[2]]), envir=dat)
  plot(fmla, ..., data=dat, type="l", xlab="", lwd=2, ylim=c(0, max(yvals)))
}


###################################################
### code chunk number 4: updates.Rnw:68-73
###################################################
getOption("SweaveHooks")[["fig"]]()
Plot((Rlines + srclines)/1000 ~ date, ylab="Lines of code (x 1000)", 
     main="Spatstat growth")
lines(srclines/1000 ~ date, data=z)
text(as.Date("2013-01-01"), 9.5, "C code")
text(as.Date("2013-01-01"), 50, "R code")


###################################################
### code chunk number 5: updates.Rnw:1390-1394
###################################################
nbugs <- nrow(news(grepl("^BUG", Category), 
                   package="spatstat"))
nbugssince <- nrow(news(Version > "1.42-0" & grepl("^BUG", Category), 
                   package="spatstat"))


###################################################
### code chunk number 6: updates.Rnw:1400-1401 (eval = FALSE)
###################################################
## news(grepl("^BUG", Category), package="spatstat")


