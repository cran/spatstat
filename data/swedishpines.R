require(spatstat, quietly=TRUE, save=FALSE)
swedishpines <- scanpp("swedishpines.tab", owin(c(0,96),c(0,100)),
                       dir=spatstat.rawdata.location())
