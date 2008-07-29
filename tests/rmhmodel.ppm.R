#
# tests of rmhmodel.ppm
#
require(spatstat)
data(cells)

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


