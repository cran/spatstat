#
#
#    pairwise.S
#
#    $Revision: 1.5 $	$Date: 2006/10/09 02:39:56 $
#
#    Pairwise()    create a user-defined pairwise interaction process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

Pairwise <- function(pot, name = "user-defined pairwise interaction process",
                     par = NULL, parnames=NULL,
                     printfun) {

  fop <- names(formals(pot))
  if(!identical(all.equal(fop, c("d", "par")), TRUE)
     && !identical(all.equal(fop, c("d", "tx", "tu", "par")), TRUE))
    stop(paste("Formal arguments of pair potential function",
               sQuote("pot"),
               "must be either (d, par) or (d, tx, tu, par)"))

  if(missing(printfun))
    printfun <- function(self) {
           cat(paste(self$name, "\n"))
           cat("Potential function:\n")
           print(self$pot)
         }

  out <- 
  list(
         name     = name,
         family   = pairwise.family,
         pot      = pot,
         par      = par,
         parnames = parnames,
         init     = NULL,
         update   = NULL,  
         print    = printfun
  )
  class(out) <- "interact"
  return(out)
}



