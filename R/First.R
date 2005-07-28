.First.lib <- function(lib, pkg) {
    library.dynam("spatstat", pkg, lib)
    v <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                  fields="Version")
    cat(paste("\nspatstat", v, "\n"))
    cat("Type \"help(spatstat)\" for information\n")
}


