#
#
#    pairwise.S
#
#    $Revision: 1.3 $	$Date: 2002/04/07 10:44:48 $
#
#    Pairwise()    create a user-defined pairwise interaction process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

Pairwise <- function(pot, name = "user-defined pairwise interaction process",
                     par = NULL, parnames=NULL,
                     printfun) {

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



