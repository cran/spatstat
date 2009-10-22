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
