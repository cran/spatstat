#
#
#    saturated.S
#
#    $Revision: 1.3 $	$Date: 2005/07/27 07:00:59 $
#
#    Saturated pairwise process with user-supplied potential
#
#    Saturated()  create a saturated pairwise process
#                 [an object of class 'interact']
#                 with user-supplied potential
#	
#
# -------------------------------------------------------------------
#	

Saturated <- function(pot, name) {
  if(missing(name))
    name <- "Saturated process with user-defined potential"
  
  fop <- names(formals(pot))
  if(!identical(all.equal(fop, c("d", "par")), TRUE)
     && !identical(all.equal(fop, c("d", "tx", "tu", "par")), TRUE))
    stop(paste("Formal arguments of pair potential function \`pot\'",
               "must be either (d, par) or (d, tx, tu, par)"))

  out <- 
  list(
         name     = name,
         family    = pairsat.family,
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
