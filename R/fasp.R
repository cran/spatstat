#
#	fasp.R
#
#	$Revision: 1.3 $	$Date: 2002/01/18 06:52:33 $
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
        included <- rep(F, length(x$fns))
        w <- as.vector(x$which[I,J])
        included[w] <- T

        # determine positions in shortened lists
        newk <- cumsum(included)

        # assign result
        Y <- list()
        Y$fns <- x$fns[included]
        Y$titles <- x$titles[included]
        xdf <- x$default.formula
        Y$default.formula <- if(any(class(xdf)=="formula")) xdf else xdf[included]
        oldwhich <- x$which[I,J,drop=F]
        newwhich <- newk[oldwhich]
        Y$which <- matrix(newwhich, ncol=ncol(oldwhich), nrow=nrow(oldwhich))
        Y$dataname <- x$dataname
        Y$title <- x$title
        class(Y) <- "fasp"
        
        return(Y)
}
