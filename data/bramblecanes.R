require(spatstat, quietly=TRUE, save=FALSE)
bramblecanes <- scanpp("bramblecanes.tab", unit.square(), multitype=TRUE,
                       dir=spatstat.rawdata.location())
