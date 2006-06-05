#
#	exactPdt.S
#	S function exactPdt() for exact distance transform of pixel image
#
#	$Revision: 4.5 $	$Date: 2006/06/01 02:29:42 $
#

"exactPdt"<-
function(w)
{
#
        verifyclass(w, "owin")
        if(w$type != "mask")
          stop("Input must be a window of type \'mask\'")
#	
	nr <- w$dim[1]
	nc <- w$dim[2]
# pad out the input image with a margin of width 1 on all sides
	x <- w$m
	x <- cbind(FALSE, x, FALSE)
	x <- rbind(FALSE, x, FALSE)
#	
	res <- .C("ps_exact_dt_R",
		as.double(w$xrange[1]),
		as.double(w$yrange[1]),
		as.double(w$xrange[2]),
		as.double(w$yrange[2]),
		nr = as.integer(nr),
		nc = as.integer(nc),
		as.logical(t(x)),
		dd = as.double (matrix(0, ncol = nc + 2, nrow = nr + 2)),
		rr = as.integer(matrix(0, ncol = nc + 2, nrow = nr + 2)),
		cc = as.integer(matrix(0, ncol = nc + 2, nrow = nr + 2)),
		bb = as.double (matrix(0, ncol = nc + 2, nrow = nr + 2)),
                PACKAGE="spatstat"
		)
	dist <- matrix(res$dd, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
	rows <- matrix(res$rr, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
	cols <- matrix(res$cc, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
	bdist<- matrix(res$bb, ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]

        # convert from C to S
        rows <- rows + 1
        cols <- cols + 1
        
	return(list(d=dist,row=rows,col=cols,b=bdist, w=w))
}
