require(spatstat, quietly=TRUE, save=FALSE)
nztrees <- scanpp("nztrees.tab", owin(c(0,153), c(0,95)),
                  dir=spatstat.rawdata.location())
