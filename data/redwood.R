require(spatstat, quietly=TRUE, save=FALSE)
redwood <- scanpp("redwood.tab", owin(c(0,1), c(-1,0)),
                  dir=spatstat.rawdata.location())
