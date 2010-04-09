# check kstest with strange data
require(spatstat)
data(ants)
# Marked point patterns with some marks not represented
AC <- split(ants, un=FALSE)$Cataglyphis
AM <- split(ants, un=FALSE)$Messor
DM <- distmap(AM)
# should produce a warning, rather than a crash:
kstest(AC, DM)
# should be OK:
kstest(unmark(AC), DM)
