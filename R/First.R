.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    cat("spatstat 1.3-1\n")
    cat("Type \"demo(spatstat)\" for a demonstration\n")
    locn <- paste(Sys.getenv("R_HOME"),"/library/spatstat/doc", sep="")
    cat(paste("See the Introduction and Quick Reference in", locn, "\n"))
}

