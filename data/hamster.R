require(spatstat, quietly=TRUE, save=FALSE)
hamster <- scanpp("hamster.tab", owin(), multitype=TRUE,
                  dir=spatstat.rawdata.location())
