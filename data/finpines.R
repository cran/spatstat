require(spatstat, quietly=TRUE, save=FALSE)
.fpath <- paste(spatstat.rawdata.location(),
                 "finpines.tab", sep=.Platform$file.sep)
.fp <- read.table(.fpath, header=TRUE)
finpines <- ppp(.fp$x, .fp$y, c(-5,5), c(-8,2), marks=.fp$height)
finpines.extra <- list(diameter=.fp$diameter)
rm(.fp, .fpath)

