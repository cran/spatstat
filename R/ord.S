#
#
#    ord.S
#
#    $Revision: 1.2 $	$Date: 2001/08/07 11:52:17 $
#
#    Ord process with user-supplied potential
#
#    Ord()  create an instance of the Ord process
#                 [an object of class 'interact']
#                 with user-supplied potential
#	
#
# -------------------------------------------------------------------
#	

Ord <- function(pot, name) {
  if(missing(name))
    name <- "Ord process with user-defined potential"
  
  out <- 
  list(
         name     = name,
         family    = ord.family,
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
