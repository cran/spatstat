#
#	tryKcross.S
#
#	Test the routine Kcross()
#
#	$Revision: 4.2 $ $Date: 2001/11/13 04:31:36 $
#
################################################################
#
"try.Kcross"<-
function(niter = 20, lambda1 = 25, lambda2 = 25, r = seq(0, 1, 0.02), R=0.2)
{
	lambda <- lambda1 + lambda2
	probs <- c(lambda1,lambda2)/lambda
	
	k <- matrix(0, nrow = niter, ncol = length(r))
	k2 <- k
	cat("computing realisation ")
	for(i in 1:niter) {
		cat(paste(i,", ", sep=""))

                X <- rpoispp(lambda)
		X$marks <- factor(sample(1:2, X$n, prob=probs, replace=T))
		out <- Kcross(X, "1", "2", r)
		k[i,  ] <- out$border
		k2[i,  ] <- out$bord.modif
	}
	cat("\n")

        # restrict to r in [0,R]
        ok <- (r <= R)
        r <- r[ok]
        k <- k[, ok]
        k2 <- k2[, ok]
        
	truek <- pi * r^2
	rrange <- range(r)
	krange <- range(c(k, k2, truek), na.rm=T)
	
	plot(rrange, krange, xlab = "r", ylab = "Kcross", 
				main = "border method", type = "n")
	for(i in 1:niter) {
		lines(r, k[i,  ])
	}
	plot(krange, krange, xlab = "true Kcross", ylab = "estimated Kcross", 
				main = "border method", type = "n")
	for(i in 1:niter) {
		lines(truek, k[i,  ])
	}
	plot(rrange, krange, xlab = "r", ylab = "Kcross", 
				main = "border (modified)", type = "n")
	for(i in 1:niter) {
		lines(r, k2[i,  ])
	}
	plot(krange, krange, xlab = "true Kcross", ylab = "estimated Kcross", 
				main = "border (modified)", type = "n")
	for(i in 1:niter) {
		lines(truek, k2[i,  ])
	}
	invisible(NULL)
}
