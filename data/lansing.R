require(spatstat, quietly=TRUE, save=FALSE)
lansing <- scanpp("lansing.tab", unit.square(), multitype=TRUE,
                  dir=spatstat.rawdata.location())
