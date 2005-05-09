#	Jest.S
#
#	Usual invocation to compute J function
#	if F and G are not required 
#
#	$Revision: 4.6 $	$Date: 2005/05/05 01:16:42 $
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

        rslt <- data.frame(r=FF$r, theo=theo, rs=Jrs, km=Jkm, un=Jun)
# convert to class "fv"
        alim <- range(rslt$r[FF$km <= 0.9])
        labl <- c("r", "Jpois(r)", "Jbord(r)", "Jkm(r)",
                  "Jun(r)")
        desc <- c("distance argument r",
                  "theoretical Poisson J(r) = 1",
                  "border corrected estimate of J(r)",
                  "Kaplan-Meier estimate of J(r)",
                  "uncorrected estimate of J(r)")
        Z <- fv(rslt, "r", "J(r)", "km", cbind(km, rs, un, 1) ~ r, alim, labl, desc)

# add more info        
        attr(Z, "F") <- FF
        attr(Z, "G") <- G

        return(Z)
}

