#	Jest.S
#
#	Usual invocation to compute J function
#	if F and G are not required 
#
#	$Revision: 4.8 $	$Date: 2005/07/15 01:56:41 $
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

        rslt <- data.frame(r=FF$r, theo=theo, un=Jun, rs=Jrs, km=Jkm)
# convert to class "fv"
        alim <- range(rslt$r[FF$km <= 0.9])
        labl <- c("r", "Jpois(r)", "Jun(r)", "Jbord(r)", "Jkm(r)")
        desc <- c("distance argument r",
                  "theoretical Poisson J(r) = 1",
                  "uncorrected estimate of J(r)",
                  "border corrected estimate of J(r)",
                  "Kaplan-Meier estimate of J(r)")
        Z <- fv(rslt, "r", substitute(J(r), NULL), "km", . ~ r, alim, labl, desc)

# add more info        
        attr(Z, "F") <- FF
        attr(Z, "G") <- G

        return(Z)
}

