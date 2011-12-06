#
# versions.R
#
# version numbers
#
# $Revision: 1.5 $  $Date: 2011/11/21 08:53:49 $
#
#####################


# Extract version string from ppm object

versionstring.ppm <- function(object) {
  verifyclass(object, "ppm")
  v <- object$version
  if(is.null(v) || !is.list(v))
    v <- list(major=1, minor=3, release=4)
  vs <- paste(v$major, ".", v$minor, "-", v$release, sep="")
  return(vs)
}

# Extract version string from interact object

versionstring.interact <- function(object) {
  verifyclass(object, "interact")
  v <- object$version
  return(v)  # NULL before 1.11-0
}

# Get version number of current spatstat installation

versionstring.spatstat <- function() {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat"),
                 fields="Version")
  return(as.character(vs))
}

# Extract major and minor versions only.

majorminorversion <- function(v) {
  vp <- package_version(v)
  return(package_version(paste(vp$major, vp$minor, sep=".")))
}

