#
#	exactPdt.S
#	S function exactPdt() for exact distance transform of pixel image
#
#	$Revision: 4.4 $	$Date: 2004/11/15 19:29:16 $
#

"exactPdt"<-
function(im)
{
#
        verifyclass(im, "owin")
        if(im$type != "mask")
          stop("Input must be a window of type \'mask\'")
#	
	nr <- im$dim[1]
	nc <- im$dim[2]
# pad out the input image with a margin of width 1 on all sides
	x <- im$m
	x <- cbind(FALSE, x, FALSE)
	x <- rbind(FALSE, x, FALSE)
#	
	res <- .C("ps_exact_dt_R",
		as.double(im$xrange[1]),
		as.double(im$yrange[1]),
		as.double(im$xrange[2]),
		as.double(im$yrange[2]),
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
        
	out <- im
	out$m <- NULL
	out <- append(out, list(d=dist,row=rows,col=cols,b=bdist))
	invisible(out)
}
