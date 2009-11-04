.First.lib <- function(lib, pkg) {
  library.dynam("spatstat", pkg, lib)
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkg),
                fields="Version")
  cat(paste("spatstat", v,
            " (Type", sQuote("help(spatstat)"), "for information)\n"))
  reset.spatstat.options()
  invisible(NULL)
}


