### R code from vignette source 'bugfixes.Rnw'

###################################################
### code chunk number 1: bugfixes.Rnw:20-24
###################################################
library(spatstat)
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: bugfixes.Rnw:38-42
###################################################
nbugs <- nrow(news(grepl("^BUG", Category), 
                   package="spatstat"))
nbugssince <- nrow(news(Version > "1.42-0" & grepl("^BUG", Category), 
                   package="spatstat"))


###################################################
### code chunk number 3: bugfixes.Rnw:58-59 (eval = FALSE)
###################################################
## bugfixes


###################################################
### code chunk number 4: bugfixes.Rnw:63-64 (eval = FALSE)
###################################################
## bugfixes(sinceversion="1.50-0")


###################################################
### code chunk number 5: bugfixes.Rnw:68-69 (eval = FALSE)
###################################################
## bugfixes(sincedate="2017-06-30")


###################################################
### code chunk number 6: bugfixes.Rnw:72-73 (eval = FALSE)
###################################################
## bugfixes("book")


###################################################
### code chunk number 7: bugfixes.Rnw:76-77 (eval = FALSE)
###################################################
## bugfixes("all")


