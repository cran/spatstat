#
#
#   allstats.R
#
#   $Revision: 1.5 $   $Date: 2001/11/06 07:09:17 $
#
#
allstats <- function(pp,dataname=NULL,verb=F) {
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

# estimate F, G and J 
  if(verb) cat("Calculating F, G, J ...")
  Jout <- Jest(pp,eps=NULL,breaks=brks)
  if(verb) cat("ok.\n")

# extract empty space function F
  Fout <- attr(Jout, "F")
  fns[[1]] <- Fout[c("r","km","rs","raw","hazard","theo")]
  titles[[1]] <- "F function"
  if(verb) cat("F done.\n")

# extract Nearest neighbour distance distribution function G
  Gout <- attr(Jout, "G")
  fns[[2]] <- Gout[c("r","km","rs","hazard","theo")]
  titles[[2]] <- "G function"
  if(verb) cat("G done.\n")

# extract J function
  fns[[3]] <- Jout[c("r","km","rs","un","theo")]
  titles[[3]] <- "J function"
  if(verb) cat("J done.\n")

# compute second moment function K
  fns[[4]] <- Kest(pp,eps=NULL,breaks=brks)[c("r","border","theo")]
  titles[[4]] <- "K function"
  if(verb) cat("K done.\n")

# wrap into 'fasp' object
  
  deform <- list(cbind(km,theo)~r,cbind(km,theo)~r,
               cbind(km,theo)~r,cbind(border,theo)~r)
  witch <- matrix(1:4,2,2,byrow=T)

  if(is.null(dataname))
    dataname <- deparse(substitute(pp))
  title <- paste("Four summary functions for ",
              	dataname,".",sep="")

  rslt <- list(fns=fns,titles=titles,default.formula=deform,which=witch,
             dataname=dataname,title=title)
  class(rslt) <- "fasp"
  rslt
}
