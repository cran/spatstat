require(spatstat, quietly=TRUE, save=FALSE)
longleaf <- scanpp("longleaf.tab", owin(c(0,200),c(0,200)),
                   dir=spatstat.rawdata.location())
