.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    cat("spatstat 1.2-1\n")
}

