require(spatstat, quietly=TRUE, save=FALSE)
simdat <- scanpp("simdat.tab", owin(c(0,10), c(0,10)),
                 dir=spatstat.rawdata.location())
