require(spatstat, quietly=TRUE, save=FALSE)
ganglia <- scanpp("ganglia.tab", owin(c(0,1), c(0,0.7533)), multitype=TRUE,
                  dir=spatstat.rawdata.location())
# relabel mark values
.mks <- ganglia$marks
levels(.mks) <- c("off", "on")
ganglia$marks <- .mks
rm(.mks)

