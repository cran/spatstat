#
#   resolve.defaults.R
#
#  $Revision: 1.1 $ $Date: 2004/08/30 04:58:34 $
#
# Resolve conflicts between several sets of defaults
# Usage:
#     resolve.defaults(list1, list2, list3, .......)
# where the earlier lists have priority 
#
resolve.defaults <- function(...) {
  arglist <- list(...)
  argue <- list()
  if((n <- length(arglist)) > 0)  {
    for(i in seq(n))
      argue <- append(argue, arglist[[i]])
  }
  if(!is.null(nam <- names(argue))) {
    named <- (nam != "")
    arg.unnamed <- argue[!named]
    arg.named <-   argue[named]
    if(any(discard <- duplicated(names(arg.named)))) 
      arg.named <- arg.named[!discard]
    argue <- append(arg.unnamed, arg.named)
  }
  return(argue)
}



  
