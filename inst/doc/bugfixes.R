### R code from vignette source 'bugfixes.Rnw'

###################################################
### code chunk number 1: bugfixes.Rnw:17-23
###################################################
library(spatstat)
x <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
              fields = c("Version", "Date"))
sversion <- as.character(x[,"Version"])
sdate    <- as.character(x[,"Date"])
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: bugfixes.Rnw:37-39
###################################################
nbugs      <- nrow(bugfixes("all",  show=FALSE))
nbugssince <- nrow(bugfixes("book", show=FALSE))


###################################################
### code chunk number 3: bugfixes.Rnw:55-56 (eval = FALSE)
###################################################
## bugfixes


###################################################
### code chunk number 4: bugfixes.Rnw:60-61 (eval = FALSE)
###################################################
## bugfixes(sinceversion="1.50-0")


###################################################
### code chunk number 5: bugfixes.Rnw:65-66 (eval = FALSE)
###################################################
## bugfixes(sincedate="2017-06-30")


###################################################
### code chunk number 6: bugfixes.Rnw:69-70 (eval = FALSE)
###################################################
## bugfixes("book")


###################################################
### code chunk number 7: bugfixes.Rnw:73-74 (eval = FALSE)
###################################################
## bugfixes("all")


###################################################
### code chunk number 8: bugfixes.Rnw:85-103
###################################################
getstuff <- function(pkg) {
  x <- read.dcf(file=system.file("DESCRIPTION", package=pkg),
                fields=c("Version", "Date"))
  xversion <- as.character(x[,"Version"])
  xdate    <- as.character(x[,"Date"])
  data.frame(date=as.Date(xdate), package=pkg, version=xversion)
}
vtable <- do.call(rbind,
                  lapply(c("spatstat.utils",
                           "spatstat.data", 
                           "spatstat.sparse",
                           "spatstat.geom",
                           "spatstat.random",
                           "spatstat.explore",
                           "spatstat.model",
                           "spatstat.linnet",
                           "spatstat"),
                         getstuff))


###################################################
### code chunk number 9: bugfixes.Rnw:109-110
###################################################
print(vtable, row.names=FALSE)


