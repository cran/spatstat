.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    v <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                  fields="Version")
    cat(paste("spatstat", v, "\n"))
    cat("Type \"demo(spatstat)\" for a demonstration\n")
    locn <- system.file("doc", package="spatstat")
    cat(paste("See the Introduction and Quick Reference in\n", locn, "\n"))
}


