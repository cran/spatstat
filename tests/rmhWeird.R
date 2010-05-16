# strange boundary cases

require(spatstat)

   if(!exists("nv"))
     nv <- 0
   if(!exists("nr"))
     nr   <- 5e3

   # Poisson process
   cat("Poisson\n")
   modP <- list(cif="poisson",par=list(beta=10), w = square(3))
   XP <- rmh(model = modP,
             start = list(n.start=25),
             control=list(nrep=nr,nverb=nv))

   # Poisson process case of Strauss
   cat("\nPoisson case of Strauss\n")
   modPS <- list(cif="strauss",par=list(beta=10,gamma=1,r=0.7), w = square(3))
   XPS <- rmh(model=modPS,
              start=list(n.start=25),
              control=list(nrep=nr,nverb=nv))
   
   # Strauss with zero intensity
   cat("\nStrauss with zero intensity\n")
   mod0S <- list(cif="strauss",par=list(beta=0,gamma=0.6,r=0.7), w = square(3))
   X0S   <- rmh(model=mod0S,start=list(n.start=80),
                     control=list(nrep=nr,nverb=nv))
   stopifnot(X0S$n == 0)

   # Poisson with zero intensity
   cat("\nPoisson with zero intensity\n")
   mod0P <- list(cif="poisson",par=list(beta=0), w = square(3))
   X0P <- rmh(model = mod0P,
             start = list(n.start=25),
             control=list(nrep=nr,nverb=nv))


   # Poisson conditioned on zero points
   cat("\nPoisson conditioned on zero points\n")
   modp <- list(cif="poisson",
                 par=list(beta=2), w = square(10))
   Xp <- rmh(modp, start=list(n.start=0), control=list(p=1, nrep=nr))
   stopifnot(Xp$n == 0)

   # Multitype Poisson conditioned on zero points
   cat("\nMultitype Poisson conditioned on zero points\n")
   modp2 <- list(cif="poisson",
                 par=list(beta=2), types=letters[1:3], w = square(10))
   Xp2 <- rmh(modp2, start=list(n.start=0), control=list(p=1, nrep=nr))
   stopifnot(is.marked(Xp2))
   stopifnot(Xp2$n == 0)

   # Multitype Poisson conditioned on zero points of each type
   cat("\nMultitype Poisson conditioned on zero points of each type\n")
   Xp2fix <- rmh(modp2, start=list(n.start=c(0,0,0)),
                 control=list(p=1, fixall=TRUE, nrep=nr))
   stopifnot(is.marked(Xp2fix))
   stopifnot(Xp2fix$n == 0)
    
