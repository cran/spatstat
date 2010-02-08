.First.lib <- function(lib, pkg) {
  library.dynam("spatstat", pkg, lib)
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkg),
                fields="Version")
  cat("spatstat", v, "\n")
  cat("Type", sQuote("help(spatstat)"), "for an overview of spatstat\n")
  cat("    ", sQuote("latest.news()"), "for news on latest version\n")
  cat("    ", sQuote("licence.polygons()"),
      "for licence information on polygon calculations\n")
  reset.spatstat.options()
  invisible(NULL)
}


