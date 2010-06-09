# Things which should cause an error
require(spatstat)

if(!exists("nv"))
  nv <- 0
if(!exists("nr"))
  nr   <- 1e3

# Strauss with zero intensity and p = 1
mod0S <- list(cif="strauss",par=list(beta=0,gamma=0.6,r=0.7), w = square(3))
out <- try(X0S   <- rmh(model=mod0S,start=list(n.start=80),
               control=list(p=1,nrep=nr,nverb=nv),verbose=FALSE))
if(!inherits(out, "try-error"))
  stop("Error not trapped (Strauss with zero intensity and p = 1) in tests/rmhErrors.R")


