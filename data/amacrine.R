amacrine <- scanpp("amacrine.tab", owin(c(0,1.6065),c(0,1)), multitype=TRUE,
                   dir=spatstat.rawdata.location())
# relabel mark values
.mks <- amacrine$marks
levels(.mks) <- c("off", "on")
amacrine$marks <- .mks
rm(.mks)

