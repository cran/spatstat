# Test examples for rmh.default
# run to reasonable length
# and with tests for validity added
# ----------------------------------------------------

require(spatstat)

if(!exists("nr"))
   nr   <- 5e3

if(!exists("nv"))
  nv <- 0
   
   # Strauss process.
   mod01 <- list(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))
   X1.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process, conditioning on n = 80:
   X2.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(p=1,nrep=nr,nverb=nv))
   stopifnot(X2.strauss$n == 80)

   # Hard core process:
   mod02 <- list(cif="hardcore",par=list(beta=2,hc=0.7),w=c(0,10,0,10))
   X3.hardcore <- rmh(model=mod02,start=list(n.start=60),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process equal to pure hardcore:
   mod02 <- list(cif="strauss",par=list(beta=2,gamma=0,r=0.7),w=c(0,10,0,10))
   X3.strauss <- rmh(model=mod02,start=list(n.start=60),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process in a polygonal window.
   x     <- c(0.55,0.68,0.75,0.58,0.39,0.37,0.19,0.26,0.42)
   y     <- c(0.20,0.27,0.68,0.99,0.80,0.61,0.45,0.28,0.33)
   mod03 <- list(cif="strauss",par=list(beta=2000,gamma=0.6,r=0.07),
                w=owin(poly=list(x=x,y=y)))
   X4.strauss <- rmh(model=mod03,start=list(n.start=90),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process in a polygonal window, conditioning on n = 42.
   X5.strauss <- rmh(model=mod03,start=list(n.start=42),
                     control=list(p=1,nrep=nr,nverb=nv))
   stopifnot(X5.strauss$n == 42)

   # Strauss process, starting off from X4.strauss, but with the
   # polygonal window replace by a rectangular one.  At the end,
   # the generated pattern is clipped to the original polygonal window.
   xxx <- X4.strauss
   xxx$window <- as.owin(c(0,1,0,1))
   X6.strauss <- rmh(model=mod03,start=list(x.start=xxx),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss with hardcore:
   mod04 <- list(cif="straush",par=list(beta=2,gamma=0.2,r=0.7,hc=0.3),
                w=c(0,10,0,10))
   X1.straush <- rmh(model=mod04,start=list(n.start=70),
                     control=list(nrep=nr,nverb=nv))
   
   # Another Strauss with hardcore (with a perhaps surprising result):
   mod05 <- list(cif="straush",par=list(beta=80,gamma=0.36,r=45,hc=2.5),
                w=c(0,250,0,250))
   X2.straush <- rmh(model=mod05,start=list(n.start=250),
                     control=list(nrep=nr,nverb=nv))
   
   # Pure hardcore (identical to X3.strauss).
   mod06 <- list(cif="straush",par=list(beta=2,gamma=1,r=1,hc=0.7),
                w=c(0,10,0,10))
   X3.straush <- rmh(model=mod06,start=list(n.start=60),
                     control=list(nrep=nr,nverb=nv))
   
   # Soft core:
   w    <- c(0,10,0,10)
   mod07 <- list(cif="sftcr",par=list(beta=0.8,sigma=0.1,kappa=0.5),
                w=c(0,10,0,10))
   X.sftcr <- rmh(model=mod07,start=list(n.start=70),
                  control=list(nrep=nr,nverb=nv))
   
   # Diggle, Gates, and Stibbard:
   mod12 <- list(cif="dgs",par=list(beta=3600,rho=0.08),w=c(0,1,0,1))
   X.dgs <- rmh(model=mod12,start=list(n.start=300),
                control=list(nrep=nr,nverb=nv))
   
   # Diggle-Gratton:
   mod13 <- list(cif="diggra",
                 par=list(beta=1800,kappa=3,delta=0.02,rho=0.04),
                 w=square(1))
   X.diggra <- rmh(model=mod13,start=list(n.start=300),
                   control=list(nrep=nr,nverb=nv))
   
   # Geyer:
   mod14 <- list(cif="geyer",par=list(beta=1.25,gamma=1.6,r=0.2,sat=4.5),
                 w=c(0,10,0,10))
   X1.geyer <- rmh(model=mod14,start=list(n.start=200),
                   control=list(nrep=nr,nverb=nv))
   
   # Geyer; same as a Strauss process with parameters
   # (beta=2.25,gamma=0.16,r=0.7):
   
   mod15 <- list(cif="geyer",par=list(beta=2.25,gamma=0.4,r=0.7,sat=10000),
                 w=c(0,10,0,10))
   X2.geyer <- rmh(model=mod15,start=list(n.start=200),
                   control=list(nrep=nr,nverb=nv))
   
   mod16 <- list(cif="geyer",par=list(beta=8.1,gamma=2.2,r=0.08,sat=3))
   data(redwood)
   X3.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(periodic=TRUE,nrep=nr,nverb=nv))
   
   # Geyer, starting from the redwood data set, simulating
   # on a torus, and conditioning on n:
   X4.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(p=1,periodic=TRUE,nrep=nr,nverb=nv))

   # Lookup (interaction function h_2 from page 76, Diggle (2003)):
      r <- seq(from=0,to=0.2,length=101)[-1] # Drop 0.
      h <- 20*(r-0.05)
      h[r<0.05] <- 0
      h[r>0.10] <- 1
      mod17 <- list(cif="lookup",par=list(beta=4000,h=h,r=r),w=c(0,1,0,1))
      X.lookup <- rmh(model=mod17,start=list(n.start=100),
                      control=list(nrep=nr,nverb=nv))
                   
   # Strauss with trend
   tr <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
   beta <- 0.3
   gmma <- 0.5
   r    <- 45
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   mod17 <- list(cif="strauss",par=list(beta=beta,gamma=gmma,r=r),w=c(0,250,0,250),
                 trend=tr3)
   X1.strauss.trend <- rmh(model=mod17,start=list(n.start=90),
                           control=list(nrep=nr,nverb=nv))

