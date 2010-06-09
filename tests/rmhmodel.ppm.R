#
# tests of rmhmodel.ppm
#
require(spatstat)
data(cells)
data(amacrine)

f <- ppm(cells)
m <- rmhmodel(f)

f <- ppm(cells, ~x)
m <- rmhmodel(f)

f <- ppm(cells, ~1, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells, ~1, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells, ~1, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells, ~1, DiggleGratton(0.05,0.1))
m <- rmhmodel(f)

f <- ppm(cells, ~1, Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

f <- ppm(cells, ~1, Geyer(0.07,2))
m <- rmhmodel(f)

f <- ppm(cells, ~1, BadGey(c(0.07,0.1,0.13),2))
m <- rmhmodel(f)

f <- ppm(cells, ~1, PairPiece(r = c(0.05, 0.1, 0.2)))
m <- rmhmodel(f)

f <- ppm(cells, ~1, AreaInter(r=0.06))
m <- rmhmodel(f)

# multitype

r <- matrix(0.07, 2, 2)
f <- ppm(amacrine, ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

h <- matrix(min(nndist(amacrine))/2, 2, 2)
f <- ppm(amacrine, ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

diag(r) <- NA
diag(h) <- NA
f <- ppm(amacrine, ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

f <- ppm(amacrine, ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

# multitype data, interaction not dependent on type

f <- ppm(amacrine, ~marks, Strauss(0.05))
m <- rmhmodel(f)

# trends

f <- ppm(cells, ~x, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells, ~y, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells, ~x+y, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells, ~polynom(x,y,2), Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

# covariates

Z <- as.im(function(x,y){ x^2+y^2 }, as.owin(cells))
f <- ppm(cells, ~z, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))

Zim <- as.im(Z, as.owin(cells))
f <- ppm(amacrine, ~z, covariates=list(z=Zim))
m <- rmhmodel(f)

Z <- as.im(function(x,y){ x^2+y }, as.owin(amacrine))
f <- ppm(amacrine, ~z + marks, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))
m <- rmhmodel(f, control=list(p=1,fixall=TRUE))

Zim <- as.im(Z, as.owin(amacrine))
f <- ppm(amacrine, ~z + marks, covariates=list(z=Zim))
m <- rmhmodel(f)
