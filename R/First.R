.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    cat("spatstat version 1.1-3\n")
}

