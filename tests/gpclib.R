#
# tests/gpclib.R
#
#  Test whether polygon geometry works with gpclib enabled
#  (only testable when gpclib is installed!)
#
#  $Revision: 1.2 $  $Date: 2011/07/28 05:20:17 $
#
require(spatstat)
if(require(gpclib)) {
  spot <- spatstat.options()
  spatstat.options(gpclib=TRUE)
  example(intersect.owin)
  example(plot.owin)
  spatstat.options(spot)
}
