#
#
#    pairwise.S
#
#    $Revision: 1.2 $	$Date: 2001/08/07 11:52:17 $
#
#    Pairwise()    create a user-defined pairwise interaction process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

Pairwise <- function(pot, name) {
  if(missing(name))
    name <- "user-defined pairwise interaction process"

  out <- 
  list(
         name     = name,
         family   = pairwise.family,
         pot      = pot,
         par      = NULL,
         parnames = NULL,
         init     = NULL,
         update   = NULL,  
         print = function(self) {
           cat(paste(self$name, "\n"))
           cat("Potential function:\n")
           print(self$pot)
           invisible()
         }
  )
  class(out) <- "interact"
  return(out)
}
