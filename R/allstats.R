#
#
#   allstats.R
#
#   $Revision: 1.10 $   $Date: 2004/01/13 06:59:57 $
#
#
allstats <- function(pp,dataname=NULL,verb=FALSE) {
#
# Function allstats --- to calculate the F, G, K, and J functions
# for an unmarked point pattern.
#
  verifyclass(pp,"ppp")
  if(is.marked(pp))
    stop("This function is applicable only to unmarked patterns.\n")

# get sensible r values
  brks <- handle.r.b.args(r=NULL, breaks=NULL, pp$window, eps=NULL)
  r    <- brks$r

# initialise  
  fns <- list()
  titles <- list()
  deform <- list()
  
# estimate F, G and J 
  if(verb) cat("Calculating F, G, J ...")
  Jout <- Jest(pp,eps=NULL,breaks=brks)
  if(verb) cat("ok.\n")

# extract empty space function F
  Fout <- attr(Jout, "F")
  fns[[1]] <- Fout
  titles[[1]] <- "F function"
  deform[[1]] <- attr(Fout, "fmla")
  if(verb) cat("F done.\n")

# extract Nearest neighbour distance distribution function G
  Gout <- attr(Jout, "G")
  fns[[2]] <- Gout
  titles[[2]] <- "G function"
  deform[[2]] <- attr(Gout, "fmla")
  if(verb) cat("G done.\n")

# extract J function
  attr(Jout, "F") <- NULL
  attr(Jout, "G") <- NULL
  fns[[3]] <- Jout
  titles[[3]] <- "J function"
  deform[[3]] <- attr(Jout, "fmla")
  if(verb) cat("J done.\n")

# compute second moment function K
  fns[[4]] <- Kout <- Kest(pp,eps=NULL,breaks=brks)
  titles[[4]] <- "K function"
  deform[[4]] <- attr(Kout, "fmla")
  if(verb) cat("K done.\n")

# wrap into 'fasp' object
  
  witch <- matrix(1:4,2,2,byrow=TRUE)

  if(is.null(dataname))
    dataname <- deparse(substitute(pp))
  title <- paste("Four summary functions for ",
              	dataname,".",sep="")

  rslt <- fasp(fns, titles, deform, witch, dataname, title)
  return(rslt)
}
