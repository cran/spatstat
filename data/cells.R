require(spatstat, quietly=TRUE, save=FALSE)
cells <- scanpp("cells.tab", unit.square(),
                dir=spatstat.rawdata.location())
