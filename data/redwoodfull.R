redwoodfull <- scanpp("redwoodfull.tab", owin(c(0,540),c(0,543)),
                      dir=spatstat.rawdata.location())
# auxiliary info
.rdiag <- list(x=c(540, 85),y = c(517, 0))
.regionI <- owin(poly=list(x=c(85,540,540),y=c(0,0,517)))
.regionII <- owin(poly=list(x=c(0,85,540,540,0),y=c(0,0,517,543,543)))
# rescale to unit square
redwoodfull <- affine(redwoodfull, mat=diag(1/c(540,543)))
redwoodfull$window <- unit.square()
.rdiag <- affinexy(.rdiag, mat=diag(1/c(540,543)))
.regionI <- affine(.regionI, mat=diag(1/c(540,543)))
.regionII <- affine(.regionII, mat=diag(1/c(540,543)))
# Ripley's subset
.regionR <- owin(c(0,0.4857),c(0,0.4857))
# 
redwoodfull.extra <-
  list(diag = .rdiag,
       regionI = .regionI,
       regionII = .regionII,
       regionR = .regionR,
       plot=function(){
         plot(redwoodfull, main="Strauss's redwood data")
         lines(redwoodfull.extra$diag)
         plot(redwoodfull.extra$regionR, add=TRUE, lty=2)
         text(0.73,0.53,"Region I")
         text(0.33,0.93,"Region II")
         text(0.12,0.47,"Ripley's subset")
       })
# clean up
rm(.rdiag, .regionI, .regionII, .regionR)
