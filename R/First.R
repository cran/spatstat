.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    cat("spatstat 1.5-3\n")
    cat("Type \"demo(spatstat)\" for a demonstration\n")
    locn <- paste(.path.package(package="spatstat"),
                  "doc", sep=.Platform$file.sep)
    cat(paste("See the Introduction and Quick Reference in\n", locn, "\n"))
}

