.First.lib <- function(lib, pkg) {
  library.dynam("spatstat", pkg, lib)
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkg),
                fields="Version")
  cat(paste("\nspatstat", v, "\n"))
  cat(paste("Type", sQuote("help(spatstat)"), "for information\n"))
}


