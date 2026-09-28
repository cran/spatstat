### R code from vignette source 'NAobjects.Rnw'

###################################################
### code chunk number 1: NAobjects.Rnw:22-23
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: NAobjects.Rnw:28-36
###################################################
library(spatstat)
requireversion(spatstat.geom, "3.5-0.003")
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sgversion <- read.dcf(file = system.file("DESCRIPTION", 
                                         package = "spatstat.geom"),
                      fields = "Version")
options(useFancyQuotes=FALSE)
set.seed(42) # for repeatability


###################################################
### code chunk number 3: NAobjects.Rnw:119-121
###################################################
X <- NAobject("ppp")
X


###################################################
### code chunk number 4: NAobjects.Rnw:137-139
###################################################
pats <- solist(cells, NAobject("ppp"), redwood)
pats


###################################################
### code chunk number 5: NAobjects.Rnw:148-150
###################################################
m <- hyperframe(X=runif(3), Y=pats)
m


###################################################
### code chunk number 6: NAobjects.Rnw:160-163
###################################################
Z <- NAobject("ppp")
is.NAobject(Z)
is.NAobject(cells)


###################################################
### code chunk number 7: NAobjects.Rnw:168-169
###################################################
inherits(Z, what="NAobject")


###################################################
### code chunk number 8: NAobjects.Rnw:181-182
###################################################
is.na(pats)


###################################################
### code chunk number 9: NAobjects.Rnw:189-192
###################################################
U <- list(cells, Z, cells)
sapply(U, is.NAobject)
sapply(U, inherits, what="NAobject")


###################################################
### code chunk number 10: NAobjects.Rnw:199-201
###################################################
h <- hyperframe(z=1:3, p=pats)
h


###################################################
### code chunk number 11: NAobjects.Rnw:208-209
###################################################
is.na(h)


###################################################
### code chunk number 12: NAobjects.Rnw:221-224
###################################################
blah <- letters[1:4]
blah[2] <- NA
blah


###################################################
### code chunk number 13: NAobjects.Rnw:230-232
###################################################
is.character(blah[2])
identical(blah[2], NA_character_)


###################################################
### code chunk number 14: NAobjects.Rnw:243-246
###################################################
Y <- rpoispp(10, nsim=3)
Y[[2]] <- NA
Y


###################################################
### code chunk number 15: NAobjects.Rnw:253-254
###################################################
solist(cells, NA, redwood)


###################################################
### code chunk number 16: NAobjects.Rnw:263-267
###################################################
g <- hyperframe(A=letters[1:3], B=rpoispp(10, nsim=3), D=runif(3))
g
g[2,2] <- NA
g


###################################################
### code chunk number 17: NAobjects.Rnw:273-275
###################################################
g[3, ] <- NA
g


###################################################
### code chunk number 18: NAobjects.Rnw:282-284
###################################################
g[,2] <- NA
g


###################################################
### code chunk number 19: NAobjects.Rnw:316-318 (eval = FALSE)
###################################################
##   X <- NAobject("ppp") 
##   K <- Kest(X)


###################################################
### code chunk number 20: NAobjects.Rnw:328-330
###################################################
  X <- NAobject("ppp")
  K <- if(is.NAobject(X)) NAobject("fv") else Kest(X)


###################################################
### code chunk number 21: NAobjects.Rnw:364-368
###################################################
A <- solapply(pats, Window)
B <- anylapply(pats, Kest)
D <- solapply(pats, Kest, demote=TRUE)
E <- anylapply(pats, npoints)


###################################################
### code chunk number 22: NAobjects.Rnw:387-391
###################################################
K <- with(m, Kest(Y))
m$G <- with(m, Gest(Y))
m$u <- with(m, clarkevans.test(Y))
with(m, u$p.value)
