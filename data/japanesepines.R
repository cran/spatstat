require(spatstat, quietly=TRUE, save=FALSE)
japanesepines <- scanpp("japanesepines.tab", square(1),
                       dir=spatstat.rawdata.location())
