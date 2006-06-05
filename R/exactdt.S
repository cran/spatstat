#
#	exactdt.S
#	S function exactdt() for exact distance transform
#
#	$Revision: 4.6 $	$Date: 2006/06/01 02:18:34 $
#

"exactdt"<-
function(X, ...)
{
	verifyclass(X, "ppp")

        w <- as.mask(X$window, ...)

#
	nr <- w$dim[1]
	nc <- w$dim[2]
#	
	res <- .C("exact_dt_R",
		as.double(X$x),
		as.double(X$y),
		as.integer(X$n),
		as.double(w$xrange[1]),
		as.double(w$yrange[1]),
		as.double(w$xrange[2]),
		as.double(w$yrange[2]),
		nr = as.integer(nr),
		nc = as.integer(nc),
		d = as.double(matrix(0, ncol = nc + 2, nrow = nr + 2)),
		i = as.integer(matrix(0, ncol = nc + 2, nrow = nr + 2)),
		b = as.double(matrix(0, ncol = nc + 2, nrow = nr + 2)),
                PACKAGE="spatstat")
	dist <- matrix(res$d, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
	inde <- matrix(res$i, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
	bdry <- matrix(res$b, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]

        inde <- inde + 1    # convert from C to S indexing
        
        return(list(d = dist, i = inde, b = bdry, w=w))
}
