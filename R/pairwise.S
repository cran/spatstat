#
#
#    pairwise.S
#
#    $Revision: 1.7 $	$Date: 2007/01/11 03:36:02 $
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
         creator  = "Pairwise",
         family   = pairwise.family,
         pot      = pot,
         par      = par,
         parnames = parnames,
         init     = NULL,
         update   = function(self, ...){
           do.call(Pairwise,
                   resolve.defaults(list(...),
                                    list(pot=self$pot, name=self$name,
                                         par=self$par, parnames=self$parnames,
                                         printfun=self$print)))
         } , 
         print    = printfun,
         version  = versionstring.spatstat()
  )
  class(out) <- "interact"
  return(out)
}



