#	Jest.S
#
#	Usual invocation to compute J function
#	if F and G are not required 
#
#	$Revision: 4.4 $	$Date: 2001/11/13 01:25:49 $
#
#
#
"Jest" <-
function(X, eps=NULL, r=NULL, breaks=NULL) {
        X <- as.ppp(X)
        brks <-  handle.r.b.args(r, breaks, X$window)$val
	FF <- Fest(X, eps, breaks=brks)
	G <- Gest(X, breaks=brks)
        ratio <- function(a, b, c) {
          result <- a/b
          result[ b == 0 ] <- c
          result
        }
	Jrs <- ratio(1-G$rs, 1-FF$rs, NA)
	Jkm <- ratio(1-G$km, 1-FF$km, NA)
	Jun <- ratio(1-G$raw, 1-FF$raw, NA)
        theo <- rep(1, length(FF$r))

        rslt <- data.frame(J=Jkm,r=FF$r,km=Jkm,rs=Jrs,un=Jun,theo=theo)
        attr(rslt, "F") <- FF
        attr(rslt, "G") <- G
        rslt
}

