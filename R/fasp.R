#
#	fasp.R
#
#	$Revision: 1.4 $	$Date: 2002/05/13 12:41:10 $
#
#
#-----------------------------------------------------------------------------
#

"[.fasp" <-
"subset.fasp" <-
  function(x, I, J, drop, ...) {

        verifyclass(x, "fasp")
        
        m <- nrow(x$which)
        n <- ncol(x$which)
        
        if(missing(I)) I <- 1:m
        if(missing(J)) J <- 1:n

        # determine index subset for lists 'fns', 'titles' etc
        included <- rep(FALSE, length(x$fns))
        w <- as.vector(x$which[I,J])
        included[w] <- TRUE

        # determine positions in shortened lists
        newk <- cumsum(included)

        # assign result
        Y <- list()
        Y$fns <- x$fns[included]
        Y$titles <- x$titles[included]
        xdf <- x$default.formula
        Y$default.formula <- if(any(class(xdf)=="formula")) xdf else xdf[included]
        oldwhich <- x$which[I,J,drop=FALSE]
        newwhich <- newk[oldwhich]
        Y$which <- matrix(newwhich, ncol=ncol(oldwhich), nrow=nrow(oldwhich))
        Y$dataname <- x$dataname
        Y$title <- x$title
        class(Y) <- "fasp"
        
        return(Y)
}
