### R code from vignette source 'replicated.Rnw'

###################################################
### code chunk number 1: replicated.Rnw:29-30
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: replicated.Rnw:35-42
###################################################
library(spatstat)
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: replicated.Rnw:180-181
###################################################
waterstriders


###################################################
### code chunk number 4: replicated.Rnw:200-201
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders, main="")


###################################################
### code chunk number 5: replicated.Rnw:208-209
###################################################
summary(waterstriders)


###################################################
### code chunk number 6: replicated.Rnw:217-218
###################################################
X <- solist(rpoispp(100), rpoispp(100), rpoispp(100))


###################################################
### code chunk number 7: replicated.Rnw:223-225
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(X)
X


###################################################
### code chunk number 8: replicated.Rnw:234-235
###################################################
A <- solist(cells, density(cells), quadratcount(cells))


###################################################
### code chunk number 9: replicated.Rnw:253-254 (eval = FALSE)
###################################################
## K <- as.anylist(lapply(waterstriders, Kest))


###################################################
### code chunk number 10: replicated.Rnw:261-262 (eval = FALSE)
###################################################
## plot(K)


###################################################
### code chunk number 11: replicated.Rnw:293-294 (eval = FALSE)
###################################################
## hyperframe(...)


###################################################
### code chunk number 12: replicated.Rnw:319-321
###################################################
H <- hyperframe(X=1:3, Y=list(sin,cos,tan))
H


###################################################
### code chunk number 13: replicated.Rnw:329-334
###################################################
G <- hyperframe(X=1:3, Y=letters[1:3], Z=factor(letters[1:3]),
                W=list(rpoispp(100),rpoispp(100), rpoispp(100)),
                U=42,
                V=rpoispp(100), stringsAsFactors=FALSE)
G


###################################################
### code chunk number 14: replicated.Rnw:364-365
###################################################
simba


###################################################
### code chunk number 15: replicated.Rnw:378-379
###################################################
pyramidal


###################################################
### code chunk number 16: replicated.Rnw:385-386
###################################################
ws <- hyperframe(Striders=waterstriders)


###################################################
### code chunk number 17: replicated.Rnw:393-395
###################################################
H$X
H$Y


###################################################
### code chunk number 18: replicated.Rnw:405-407
###################################################
H$U <- letters[1:3]
H


###################################################
### code chunk number 19: replicated.Rnw:412-416
###################################################
G <- hyperframe()
G$X <- waterstriders
G$Y <- 1:3
G


###################################################
### code chunk number 20: replicated.Rnw:424-428
###################################################
H[,1]
H[2,]
H[2:3, ]
H[1,1]


###################################################
### code chunk number 21: replicated.Rnw:434-437
###################################################
H[,1,drop=TRUE]
H[1,1,drop=TRUE]
H[1,2,drop=TRUE]


###################################################
### code chunk number 22: replicated.Rnw:450-458 (eval = FALSE)
###################################################
## plot.solist(x, ..., main, arrange = TRUE, nrows = NULL, ncols = NULL,
##             main.panel = NULL, mar.panel = c(2, 1, 1, 2), hsep = 0, vsep = 0, 
##             panel.begin = NULL, panel.end = NULL, panel.args = NULL, 
##             panel.begin.args = NULL, panel.end.args = NULL, panel.vpad = 0.2, 
##             plotcommand = "plot", do.plot = TRUE, adorn.left = NULL, 
##             adorn.right = NULL, adorn.top = NULL, adorn.bottom = NULL, 
##             adorn.size = 0.2, adorn.args = list(), equal.scales = FALSE, 
##             halign = FALSE, valign = FALSE)


###################################################
### code chunk number 23: replicated.Rnw:473-474
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders, pch=16, nrows=1)


###################################################
### code chunk number 24: replicated.Rnw:493-494
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simba)


###################################################
### code chunk number 25: replicated.Rnw:506-508
###################################################
getOption("SweaveHooks")[["fig"]]()
H <- hyperframe(X=1:3, Y=list(sin,cos,tan))
plot(H$Y)


###################################################
### code chunk number 26: replicated.Rnw:521-522 (eval = FALSE)
###################################################
## plot(h, e)


###################################################
### code chunk number 27: replicated.Rnw:531-532
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }))


###################################################
### code chunk number 28: replicated.Rnw:544-546
###################################################
getOption("SweaveHooks")[["fig"]]()
H <- hyperframe(Bugs=waterstriders)
plot(H, quote(plot(Kest(Bugs), legend=FALSE)), marsize=2)


###################################################
### code chunk number 29: replicated.Rnw:559-561
###################################################
df <- data.frame(A=1:10, B=10:1)
with(df, A-B)


###################################################
### code chunk number 30: replicated.Rnw:574-575 (eval = FALSE)
###################################################
## with(h,e)


###################################################
### code chunk number 31: replicated.Rnw:585-588
###################################################
H <- hyperframe(Bugs=waterstriders)
with(H, npoints(Bugs))
with(H, distmap(Bugs))


###################################################
### code chunk number 32: replicated.Rnw:611-612
###################################################
with(simba, npoints(Points))


###################################################
### code chunk number 33: replicated.Rnw:619-621
###################################################
H <- hyperframe(Bugs=waterstriders)
K <- with(H, Kest(Bugs))


###################################################
### code chunk number 34: replicated.Rnw:629-630
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(K)


###################################################
### code chunk number 35: replicated.Rnw:635-637
###################################################
H <- hyperframe(Bugs=waterstriders)
with(H, nndist(Bugs))


###################################################
### code chunk number 36: replicated.Rnw:643-644
###################################################
with(H, min(nndist(Bugs)))


###################################################
### code chunk number 37: replicated.Rnw:656-657
###################################################
simba$Dist <- with(simba, distmap(Points))


###################################################
### code chunk number 38: replicated.Rnw:670-674
###################################################
getOption("SweaveHooks")[["fig"]]()
lambda <- rexp(6, rate=1/50)
H <- hyperframe(lambda=lambda)
H$Points <- with(H, rpoispp(lambda))
plot(H, quote(plot(Points, main=paste("lambda=", signif(lambda, 4)))))


###################################################
### code chunk number 39: replicated.Rnw:680-681
###################################################
H$X <- with(H, rpoispp(50))


###################################################
### code chunk number 40: replicated.Rnw:713-715
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simba, quote(plot(density(Points), main="")), 
     main="", nrows=2, marsize=1.5)


###################################################
### code chunk number 41: replicated.Rnw:740-742
###################################################
getOption("SweaveHooks")[["fig"]]()
rhos <- with(demohyper, rhohat(Points, Image))
plot(rhos, legend=FALSE, main="")


###################################################
### code chunk number 42: replicated.Rnw:759-760 (eval = FALSE)
###################################################
## mppm(formula, data, interaction, ...)


###################################################
### code chunk number 43: replicated.Rnw:770-771 (eval = FALSE)
###################################################
## mppm(Points ~ group, simba, Poisson())


###################################################
### code chunk number 44: replicated.Rnw:804-805
###################################################
mppm(Points ~ 1, simba)


###################################################
### code chunk number 45: replicated.Rnw:812-813
###################################################
mppm(Points ~ group, simba)


###################################################
### code chunk number 46: replicated.Rnw:819-820
###################################################
mppm(Points ~ id, simba)


###################################################
### code chunk number 47: replicated.Rnw:830-831
###################################################
mppm(Points ~ Image, data=demohyper)


###################################################
### code chunk number 48: replicated.Rnw:849-850 (eval = FALSE)
###################################################
## mppm(Points ~ offset(log(Image)), data=demohyper)


###################################################
### code chunk number 49: replicated.Rnw:862-863 (eval = FALSE)
###################################################
## mppm(Points ~ log(Image), data=demop)


###################################################
### code chunk number 50: replicated.Rnw:880-881 (eval = FALSE)
###################################################
## mppm(formula, data, interaction, ..., iformula=NULL)


###################################################
### code chunk number 51: replicated.Rnw:931-932
###################################################
radii <- with(simba, mean(nndist(Points)))


###################################################
### code chunk number 52: replicated.Rnw:939-941
###################################################
Rad <- hyperframe(R=radii)
Str <- with(Rad, Strauss(R))


###################################################
### code chunk number 53: replicated.Rnw:946-948
###################################################
Int <- hyperframe(str=Str)
mppm(Points ~ 1, simba, interaction=Int)


###################################################
### code chunk number 54: replicated.Rnw:975-978
###################################################
h <- hyperframe(Y=waterstriders)
g <- hyperframe(po=Poisson(), str4 = Strauss(4), str7= Strauss(7))
mppm(Y ~ 1, data=h, interaction=g, iformula=~str4)


###################################################
### code chunk number 55: replicated.Rnw:989-990
###################################################
fit <- mppm(Points ~ 1, simba, Strauss(0.07), iformula = ~Interaction*group)


###################################################
### code chunk number 56: replicated.Rnw:1008-1009
###################################################
fit


###################################################
### code chunk number 57: replicated.Rnw:1012-1014
###################################################
co <- coef(fit)
si <- function(x) { signif(x, 4) }


###################################################
### code chunk number 58: replicated.Rnw:1025-1026
###################################################
coef(fit)


###################################################
### code chunk number 59: replicated.Rnw:1083-1084 (eval = FALSE)
###################################################
## interaction=hyperframe(po=Poisson(), str=Strauss(0.07))


###################################################
### code chunk number 60: replicated.Rnw:1089-1090 (eval = FALSE)
###################################################
## iformula=~ifelse(group=="control", po, str)


###################################################
### code chunk number 61: replicated.Rnw:1100-1101 (eval = FALSE)
###################################################
## iformula=~I((group=="control")*po) + I((group=="treatment") * str)


###################################################
### code chunk number 62: replicated.Rnw:1111-1116
###################################################
g <- hyperframe(po=Poisson(), str=Strauss(0.07))
fit2 <- mppm(Points ~ 1, simba, g, 
             iformula=~I((group=="control")*po) 
                     + I((group=="treatment") * str))
fit2


###################################################
### code chunk number 63: replicated.Rnw:1139-1142
###################################################
H <- hyperframe(W=waterstriders)
fit <- mppm(W ~ 1, H)
subfits(fit)


###################################################
### code chunk number 64: replicated.Rnw:1163-1164 (eval = FALSE)
###################################################
## subfits <- subfits.new


###################################################
### code chunk number 65: replicated.Rnw:1176-1178
###################################################
H <- hyperframe(W=waterstriders)
with(H, ppm(W))


###################################################
### code chunk number 66: replicated.Rnw:1201-1203
###################################################
fit <- mppm(P ~ x, hyperframe(P=waterstriders))
res <- residuals(fit)


###################################################
### code chunk number 67: replicated.Rnw:1213-1214
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(res)


###################################################
### code chunk number 68: replicated.Rnw:1219-1221
###################################################
getOption("SweaveHooks")[["fig"]]()
smor <- with(hyperframe(res=res), Smooth(res, sigma=4))
plot(smor, main="")


###################################################
### code chunk number 69: replicated.Rnw:1233-1236
###################################################
fit <- mppm(P ~ x, hyperframe(P=waterstriders))
res <- residuals(fit)
totres <- sapply(res, integral.msr)


###################################################
### code chunk number 70: replicated.Rnw:1242-1250
###################################################
getOption("SweaveHooks")[["fig"]]()
fit <- mppm(Points~Image, data=demohyper)
resids <- residuals(fit, type="Pearson")
totres <- sapply(resids, integral.msr)
areas <- with(demohyper, area.owin(as.owin(Points)))
df <- as.data.frame(demohyper[, "Group"])
df$resids <- totres/areas
par(mar=rep(2, 4))
plot(resids~Group, df)


###################################################
### code chunk number 71: replicated.Rnw:1274-1277
###################################################
getOption("SweaveHooks")[["fig"]]()
fit <- mppm(P ~ 1, hyperframe(P=waterstriders))
sub <- hyperframe(Model=subfits(fit))
plot(sub, quote(diagnose.ppm(Model)), main="", marsize=1.5)


###################################################
### code chunk number 72: replicated.Rnw:1290-1298
###################################################
H <- hyperframe(P = waterstriders)
fitall <- mppm(P ~ 1, H)
together <- subfits(fitall)
separate <- with(H, ppm(P))
Fits <- hyperframe(Together=together, Separate=separate)
dr <- with(Fits, unlist(coef(Separate)) - unlist(coef(Together)))
dr
exp(dr)


###################################################
### code chunk number 73: replicated.Rnw:1315-1324
###################################################
H <- hyperframe(X=waterstriders)

# Poisson with constant intensity for all patterns
fit1 <- mppm(X~1, H)
quadrat.test(fit1, nx=2)

# uniform Poisson with different intensity for each pattern
fit2 <- mppm(X ~ id, H)
quadrat.test(fit2, nx=2)


###################################################
### code chunk number 74: replicated.Rnw:1353-1354 (eval = FALSE)
###################################################
## kstest.mppm(model, covariate)
