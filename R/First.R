.First.lib <- function(lib, pkg) {
  library.dynam("spatstat", pkg, lib)
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkg),
                fields="Version")
  cat(paste("spatstat", v,
            " (Type", sQuote("help(spatstat)"), "for information)\n"))
  reset.spatstat.options()
  if(!spatstat.options("gpclib")) 
    cat(paste("\n",
              "\tNote: polygon geometry computations in spatstat\n",
              "\tdepend on the package gpclib, which has a\n",
              "\trestricted licence. It is disabled by default;\n",
              "\tto enable gpclib, type",
              sQuote("spatstat.options(gpclib=TRUE)"),
              "\n"))
  invisible(NULL)
}


