require(spatstat, quietly=TRUE, save=FALSE)
.bpath <- paste(spatstat.rawdata.location(),
                 "betacells.tab", sep=.Platform$file.sep)
.beta <- read.table(.bpath, header=TRUE)
betacells <- ppp(.beta$x, .beta$y, c(28.08, 778.08), c(16.2, 1007.02),
                marks=.beta$type)
betacells.extra <- list(area=.beta$area)
rm(.bpath, .beta)

