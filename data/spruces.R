require(spatstat, quietly=TRUE, save=FALSE)
spruces <- scanpp("spruces.tab",
                  owin(c(0,56), c(0,38)),
                  multitype=FALSE,
                  dir=spatstat.rawdata.location())



