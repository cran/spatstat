#
#   When the package is installed, this tells us 
#   the directory where the .tab files are stored
#
#   Typically data/murgatroyd.R reads data-raw/murgatroyd.tab
#   and applies special processing
#
spatstat.rawdata.location <- function() {
    locn <- paste(.path.package(package="spatstat"),
                  "data-raw", sep=.Platform$file.sep)
    return(locn)
}
