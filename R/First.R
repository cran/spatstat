.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    cat("spatstat 1.5-5\n")
    cat("Type \"demo(spatstat)\" for a demonstration\n")
    locn <- system.file("doc", package="spatstat")
    cat(paste("See the Introduction and Quick Reference in\n", locn, "\n"))
}

