#
# pppmatch.R
#
# $Revision: 1.3 $  $Date: 2006/06/30 10:05:55 $
#
# From original code by Dominic Schuhmacher
#
#
# -----------------------------------------------------------------
# The standard functions for the new class pppmatching
#
# Objects of class pppmatching consist of two point patterns pp1 and pp2,
# an adjacency matrix ((i,j)-th entry 1 if i-th point of pp1 and j-th
# point of pp2 are matched, 0 otherwise), and a string for the type
# of the matching (e.g. Wasserstein, Prohorov).
# -----------------------------------------------------------------

pppmatching <- function(X, Y, am, ty = "generic") {
   verifyclass(X, "ppp")
   verifyclass(Y, "ppp")
   am <- as.matrix(am)
   am <- apply(am, c(1,2), as.logical)
   if (dim(am)[1] != X$n || dim(am)[2] != Y$n)
      stop("Adjacency matrix does not have the right dimensions")
     
   res <- list("pp1"=X, "pp2"=Y, "am"=am, "type"=ty)
   class(res) <- "pppmatching"
   res
}

plot.pppmatching <- function(x, addmatch = NULL, main=NULL, ...) {
   if (is.null(main))
        main <- deparse(substitute(x))
   pp1 <- x$pp1
   pp2 <- x$pp2
   plot.owin(pp1$window, main = main, ...)
   here <- which(x$am, arr.ind = TRUE)
   if (!is.null(addmatch)) {
      addhere <- which(addmatch, arr.ind = TRUE)
      seg <- as.psp(from=pp1[addhere[,1]], to=pp2[addhere[,2]])
      plot(seg, add=TRUE, lty = 2, col="gray70")
   }
   if (length(here) > 0) {
     seg <- as.psp(from=pp1[here[,1]], to=pp2[here[,2]])
     plot(seg, add=TRUE, ...)
   }
   oldpar <- par("pch", "col")
   par(pch=16, col=2)
   points(x$pp1, ...)
   par(col=4)
   points(x$pp2, ...)
   par(oldpar)
   return(invisible(NULL))
}

print.pppmatching <- function(x, ...) {
   cat(x$type, "matching of two planar point patterns \n")
   cat("pp1:", x$pp1$n, "points \n")
   cat("pp2:", x$pp2$n, "points \n")
   print.owin(x$pp1$window)
   cat("pairing:", sum(x$am), "lines \n")
   return(invisible(NULL))
}

summary.pppmatching <- function(object, ...) {
   cat(paste(object$type, "matching of two planar point patterns \n"))
   cat(paste("pp1:", object$pp1$n, "points \n"))
   cat(paste("pp2:", object$pp2$n, "points \n"))
   print.owin(object$pp1$window)
   if ((npair <- sum(object$am)) == 0) cat("pairing: empty \n") 
   else {
     cat(paste("pairing:", npair, ngettext(npair, "line", "lines"), "\n \n"))
     rowsum <- apply(object$am, 1, "sum")
     colsum <- apply(object$am, 1, "sum")
     lt <- ifelse(min(rowsum) >= 1, TRUE, FALSE)
     ru <- ifelse(max(rowsum) <= 1, TRUE, FALSE)
     rt <- ifelse(min(colsum) >= 1, TRUE, FALSE)
     lu <- ifelse(max(colsum) <= 1, TRUE, FALSE)
     if (lt && ru && rt && lu)
       cat("pairing is 1-1 \n")
     else if (!lt && !ru && !rt && !lu)
       cat("pairing has no special properties \n")
     else 
       cat(paste("pairing is",
                 ifelse(lt, "left-total", ""),
                 ifelse(ru, "right-unique", ""),
                 ifelse(rt, "right-total", ""),
                 ifelse(lu, "left-unique", ""),
                 "\n") )
   }
   return(invisible(NULL))
 }


# -----------------------------------------------------------------
# Implementation of the primal-dual assignment algorithm
# -----------------------------------------------------------------
#
# pppdist is the main function (maxflow is used as a "lemma")
#
# Its arguments are 
#
# x and y of class ppp (the two point patterns for which we want to do the
#   assignment)
# the order q of the Wasserstein metric used (so q=1 for the usual metric
#   (minimization of average pairing distance), other q for its l_q generalizations,
#    l=Inf is not yet integrated, but there is a separate function ppprohorov below)
# precision and belowone have something to do with the numerical precision and
#   are somewhat experimental (but only of concern for bigger q, as long as x and
#   y are point processes on the unit square)
# show.rprimal=TRUE shows at each stage of the algorithm what the current restricted
#   primal problem and its solution are (algorithm jumps between restricted primal
#   and dual problem until the solution to the restricted primal (a partial
#   pairing of the point patterns) is a full pairing)
# timelag gives the number of seconds of pause added each time a solution to
#   the current restricted primal is found (has only an effect if show.primal=TRUE) 
#   
# Note that if show.rprimal=TRUE ("show mode"), my old R code is used (which is very
#   slow), whereas if show.rprimal=FALSE ("fast computation mode") the C program is
#   called which should solve problems for about hundred points per pattern in
#   instant.
#
# -----------------------------------------------------------------

pppdist <- function(X, Y, q=1, precision=7,
                   show.rprimal = FALSE, belowone = TRUE, timelag = 0) {
  verifyclass(X, "ppp")
  verifyclass(Y, "ppp")
  if(is.infinite(q))
    return(pppdist.prohorov(X, Y, precision))
  if (X$n != Y$n) {
    message("Total numbers of points do not match. Distance is maximal.")
    return()
  }
  if (X$n == 0) {
    message("Point patterns are empty. Distance is zero")
    return()
  }
  n <- X$n
  dfix <- crossdist(X, Y)
  fact = ifelse(belowone, sqrt(2), 1)
  d <- round((dfix/fact)^q*10^precision)

  if (show.rprimal) {
      plot(pppmatching(X, Y, matrix(FALSE, n, n)))
      # initialization of dual variables
      u <- apply(d, 1, min)
      d <- d - u
      v <- apply(d, 2, min)
      d <- d - rep(v, each=n)
      # the main loop
      feasible <- FALSE
      while (!feasible) {
         rpsol <- maxflow(d)  # rpsol = restricted primal, solution
         am <- matrix(FALSE, n, n)
         for (i in 1:n) {
            if (rpsol$assignment[i] > -1) am[i, rpsol$assignment[i]] <- TRUE
         }
         Sys.sleep(timelag)
         channelmat = (d == 0 & !am)
         plot(pppmatching(X, Y, am), addmatch = channelmat)
         # if the solution of the restricted primal is not feasible for  
         # the original primal, update dual variables
         if (min(rpsol$assignment) == -1) {
            w1 <- which(rpsol$fi_rowlab > -1)
            w2 <- which(rpsol$fi_collab == -1)
            subtractor <- min(d[w1, w2])
            d[w1,] <- d[w1,] - subtractor
            d[,-w2] <- d[,-w2] + subtractor 
         }
         # otherwise break the loop
         else {
            feasible <- TRUE
         }   
      }
      assig <- rpsol$assignment
   }
   else
   {
      res <- .C("done_R",
                as.integer(d),
                as.integer(n),
                assignment = as.integer(rep(-1,n)),
                PACKAGE="spatstat")
      assig <- res$assignment
      am <- matrix(FALSE, n, n)
      am[cbind(1:n, assig[1:n])] <- TRUE
      resdist = mean((dfix[am])^q)^(1/q)
      plot(pppmatching(X, Y, am), main = paste("p = ", q, ",  distance = ",
         format(resdist, digits=4), sep=""))
   }
   resdist = mean((dfix[am])^q)^(1/q)
   print(resdist)
   return(pppmatching(X, Y, am, paste("Wasserstein", q, sep="")))   
}   

#  
# Solution of restricted primal
# 

maxflow <- function(costm) {

  stopifnot(is.matrix(costm))
  stopifnot(nrow(costm) == ncol(costm))
  if(!all(apply(costm == 0, 1, any)))
    stop("Each row of the cost matrix must contain a zero")
  
  m <- dim(costm)[1]   # cost matrix is square m * m
  assignment <- rep(-1, m)   # -1 means no pp2-point assigned to i-th pp1-point
   # initial assignment or rowlabel <- source label (= 0) where not possible
   for (i in 1:m) {
      j <- match(0, costm[i,])
      if (!(j %in% assignment))
         assignment[i] <- j
   }
   newlabelfound <- TRUE
   while (newlabelfound) {
     rowlab <- rep(-1, m)   # -1 means no label given, 0 stands for source label
     collab <- rep(-1, m)
     rowlab <- ifelse(assignment == -1, 0, rowlab)
     # column and row labeling procedure until either breakthrough occurs
     # (which means that there is a better point assignment, i.e. one that
     # creates more point pairs than the current one (flow can be increased))
     # or no more labeling is possible
     breakthrough <- -1
     while (newlabelfound && breakthrough == -1) { 
         newlabelfound <- FALSE
         for (i in 1:m) {
            if (rowlab[i] != -1) {
               for (j in 1:m) {
                  if (costm[i,j] == 0 && collab[j] == -1) {
                     collab[j] <- i
                     newlabelfound <- TRUE
                     if (!(j %in% assignment) && breakthrough == -1)
                        breakthrough <- j
                  }
               }
            }
         }
         for (j in 1:m) {
            if (collab[j] != -1) {
               for (i in 1:m) {
                  if (assignment[i] == j && rowlab[i] == -1) {
                     rowlab[i] <- j
                     newlabelfound <- TRUE
                  }
               }
            }
         }
      }
      # if the while-loop was left due to breakthrough,
      # reassign points (i.e. redirect flow) and restart labeling procedure
      if (breakthrough != -1) {
         l <- breakthrough
         while (l != 0) {
            k <- collab[l]
            assignment[k] <- l
            l <- rowlab[k] 
         }
      }
   }
   # the outermost while-loop is left, no more labels can be given; hence
   # the maximal number of points are paired given the current restriction
   # (flow is maximal given the current graph)
   return(list("assignment"=assignment, "fi_rowlab"=rowlab, "fi_collab"=collab))  
}







# -----------------------------------------------------------------
# Implementation of brute force Prohorov distance computation
# ("the case q=Inf in ppdist")
# -----------------------------------------------------------------

pppdist.prohorov <- function(X, Y, precision=7) {
   if (X$n != Y$n) {
      message("Total numbers of points do not match. Distance is maximal.")
      return()
   }
   if (X$n == 0) {
      message("Point patterns are empty. Distance is zero")
      return()
   }
   n <- X$n
   dfix <- crossdist(X, Y)
   d <- round(dfix*10^precision)
   res <- .C("dinfty_R",
             as.integer(d),
             as.integer(n),
             assignment = as.integer(rep(-1,n)),
             PACKAGE="spatstat")
   assig <- res$assignment
   am <- matrix( FALSE, n, n)
   for (i in 1:n) am[i, assig[i]] <- TRUE
   resdist <- max(dfix[am])
   print(resdist)
   pm <- pppmatching(X, Y, am, paste("Prohorov"))
   plot(pm)
   return(pm)
}   





