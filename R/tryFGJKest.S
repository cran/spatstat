#
#	tryFGJKest.S
#
#	Test the F, G, J and K estimation routines
#
#	$Revision: 4.3 $ $Date: 2002/05/13 12:41:10 $
#
#   Note: this is not the most efficient way to calculate F, G, J and K
#         if all four are required. See the source for 'allstats'
#
################################################################
#
"try.FGJKest"<-
function(niter = 20, lambda = 25, r = seq(0, sqrt(2), 0.02), eps=0.01, slow=FALSE)
{
	Frs <- 
	Fkm <- 
	Grs <- 
	Gkm <- 
	Jrs <- 
	Jkm <- 
	Kbord <- matrix(0, nrow = niter, ncol = length(r))

	cat("computing realisation ")
	for(i in 1:niter) {
		cat(paste(i))

		pp <- rpoispp(lambda)
                cat("(")
                cat("F")
                FF <- Fest(pp, eps, r)
		Fkm[i,  ] <- FF$km
		Frs[i,  ] <- FF$rs
                cat("G")
                G <- Gest(pp, r)
		Gkm[i,  ] <- G$km
		Grs[i,  ] <- G$rs
                cat("J")
                J <- Jest(pp, eps, r)
		Jkm[i,  ] <- J$km
		Jrs[i,  ] <- J$rs
                cat("K") ;
                K <- Kest(pp, r, slow=slow)
		Kbord[i,  ] <- K$border
                cat("), ")
	}
	cat("Done.\n")

        oldpar <- par(ask=TRUE)
	plotteststuff <- function(r, mat, correction, symb, trueval) {
          	rrange <- range(r)
                mrange <- range(c(mat, trueval), na.rm=TRUE)
                plot(rrange, mrange, xlab="r", ylab=symb, sub=correction,
                     main=paste("estimated", symb), type="n")
                niter <- nrow(mat)
                for(i in 1:niter) {
                  lines(r, mat[i,  ])
                }
                lines(r, trueval, lty=2)
                dev <- mat - t(matrix(trueval, ncol=nrow(mat), nrow=ncol(mat)))
                derange <- range(c(dev,0), na.rm=TRUE)
                plot(rrange, derange, xlab="r", ylab="Deviation",
                     sub=correction,
                     main=paste("deviation of estimated", symb), type="n")
                for(i in 1:niter) {
                  lines(r, dev[i,  ])
                }
                abline(0,0,lty=2)
                bias <- apply(mat, 2, mean, na.rm=TRUE) - trueval
                brange <- range(c(bias, 0), na.rm=TRUE)
                plot(r, bias, xlab="r", ylab="Bias",
                     sub=correction, ylim=brange,
                     main=paste("Bias of estimated", symb), type="l")
                abline(0,0,lty=2)
                
                varnaok <- function(x) {
                  nbg <- is.na(x)
                  if(all(nbg))
                    NA
                  else
                    var(x[!nbg])   # should work in all dialects
                }
                sd <- sqrt(apply(mat, 2, varnaok))
# Splus 5.1:    sd <- sqrt(apply(mat, 2, var, na.method="available"))
# R:            sd <- sqrt(apply(mat, 2, var, na.rm=TRUE))
                
                plot(r, sd, xlab="r", ylab="SD",
                     sub=correction,
                     main=paste("Standard deviation of estimated", symb),
                     type="l")
                invisible(NULL)
              }

        trueK  <- pi * r^2
        trueFG <- 1 - exp( - lambda * pi * r^2)
        trueJ <- rep(1, length(r))

        plotteststuff(r, Fkm, "Kaplan-Meier", "F", trueFG)
        plotteststuff(r, Frs, "reduced sample", "F", trueFG)
        plotteststuff(r, Gkm, "Kaplan-Meier", "G", trueFG)
        plotteststuff(r, Grs, "reduced sample", "G", trueFG)
        plotteststuff(r, Jkm, "Kaplan-Meier", "J", trueJ)
        plotteststuff(r, Jrs, "reduced sample", "J", trueJ)
        plotteststuff(r, Kbord, "border method", "K", trueK)

        par(oldpar)
	invisible(NULL)
}
