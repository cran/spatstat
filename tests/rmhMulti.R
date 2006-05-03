#
#  tests of rmh, running multitype point processes
#
require(spatstat)

if(!exists("nr"))
   nr   <- 1e5

if(!exists("nv"))
   nv   <- 0

   # Multitype Poisson
   modp2 <- list(cif="poisson",
                 par=list(beta=2), types=letters[1:3], w = square(10))
   Xp2 <- rmh(modp2, start=list(n.start=0), control=list(p=1))
    
   # Multitype Strauss:
   beta <- c(0.027,0.008)
   gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
   r    <- matrix(c(45,45,45,45),2,2)
   mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gmma,radii=r),
                w=c(0,250,0,250))
   X1.straussm <- rmh(model=mod08,start=list(n.start=80),
                      control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   

   # Multitype Strauss conditioning upon the total number
   # of points being 80:
   X2.straussm <- rmh(model=mod08,start=list(n.start=80),
                      control=list(p=1,ptypes=c(0.75,0.25),nrep=nr,
                                   nverb=nv))
   stopifnot(X2.straussm$n == 80)

   # Conditioning upon the number of points of type 1 being 60
   # and the number of points of type 2 being 20:
   X3.straussm <- rmh(model=mod08,start=list(n.start=c(60,20)),
                      control=list(fixall=TRUE,p=1,ptypes=c(0.75,0.25),
                                   nrep=nr,nverb=nv))
   stopifnot(all(table(X3.straussm$marks) == c(60,20)))

   # Multitype Strauss hardcore:
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   mod09 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                iradii=r,hradii=rhc),w=c(0,250,0,250))
   X.straushm <- rmh(model=mod09,start=list(n.start=80),
                     control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))

   # Multitype Strauss hardcore with trends for each type:
   beta  <- c(0.27,0.08)
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   tr4   <- function(x,y){x <- x/250; y <- y/250;
                         exp(-0.6*x+0.5*y)}
                        # log linear trend
   mod10 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=c(0,250,0,250),
                 trend=list(tr3,tr4))
   X1.straushm.trend <- rmh(model=mod10,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),
                            nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore with trends for each type, given as images:
   bigwin <- square(250)
   i1 <- as.im(tr3, bigwin)
   i2 <- as.im(tr4, bigwin)
   mod11 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=bigwin,
                 trend=list(i1,i2))
   X2.straushm.trend <- rmh(model=mod11,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),expand=1,
                            nrep=nr,nverb=nv))


#######################################################################
############  checks on distribution of output  #######################
#######################################################################

checkp <- function(p, context, testname, failmessage, pcrit=0.01) {
  if(missing(failmessage))
    failmessage <- paste("output failed", testname)
  if(p < pcrit)
    warning(paste(context, ",",  failmessage), call.=FALSE)
  cat(paste("\n", context, ",", testname, "has p-value", signif(p,4), "\n"))
}

# Multitype Strauss code; output is multitype Poisson

beta  <- 100 * c(1,1)
ri    <- matrix(0.07, 2, 2)
gmma  <- matrix(1, 2, 2)  # no interaction
tr1   <- function(x,y){ rep(1, length(x)) }
tr2   <- function(x,y){ rep(2, length(x)) }
mod <- rmhmodel(cif="straussm",
                  par=list(beta=beta,gamma=gmma,radii=ri),
                  w=owin(),
                  trend=list(tr1,tr2))

X <- rmh(mod, start=list(n.start=0), control=list(nrep=1e6))

# The model is Poisson with intensity 100 for type 1 and 200 for type 2.
# Total number of points is Poisson (300)
# Marks are i.i.d. with P(type 1) = 1/3, P(type 2) = 2/3.

# Test whether the total intensity looks right
#
p <- ppois(X$n, 300)
p.val <- 2 * min(p, 1-p)
checkp(p.val, 
       "In multitype Poisson simulation",
       "test whether total number of points has required mean value")

# Test whether the mark distribution looks right
ta <- table(X$marks)
cat("Frequencies of marks:")
print(ta)
checkp(chisq.test(ta, p = c(1,2)/3)$p.value,
       "In multitype Poisson simulation",
       "chi-squared goodness-of-fit test for mark distribution (1/3, 2/3)")

#####
####  multitype Strauss code; fixall=TRUE;
####  output is multinomial process with nonuniform locations
####

the.context <- "In nonuniform multinomial simulation"

beta  <- 100 * c(1,1)
ri    <- matrix(0.07, 2, 2)
gmma  <- matrix(1, 2, 2)  # no interaction
tr1   <- function(x,y){ ifelse(x < 0.5, 0, 2) } 
tr2   <- function(x,y){ ifelse(y < 0.5, 1, 3) }
# cdf of these distributions
Fx1 <- function(x) { ifelse(x < 0.5, 0, ifelse(x < 1, 2 * x - 1, 1)) }
Fy2 <- function(y) { ifelse(y < 0, 0,
                           ifelse(y < 0.5, y/2,
                                  ifelse(y < 1, (1/2 + 3 * (y-1/2))/2, 1))) }
                                                               

mod <- rmhmodel(cif="straussm",
                  par=list(beta=beta,gamma=gmma,radii=ri),
                  w=owin(),
                  trend=list(tr1,tr2))

X <- rmh(mod, start=list(n.start=c(50,50)),
           control=list(nrep=1e6, expand=1, p=1, fixall=TRUE))

# The model is Poisson 
# Mean number of type 1 points = 100
# Mean number of type 2 points = 200
# Total intensity = 300
# Marks are i.i.d. with P(type 1) = 1/3, P(type 2) = 2/3

# Test whether the coordinates look OK
Y <- split(X)
X1 <- Y[[names(Y)[1]]]
X2 <- Y[[names(Y)[2]]]
checkp(ks.test(X1$y, "punif")$p.value,
       the.context,
       "Kolmogorov-Smirnov test of uniformity of y coordinates of type 1 points")
if(any(X1$x < 0.5)) {
  warning(paste(the.context, ",", 
                "x-coordinates of type 1 points are IMPOSSIBLE"), call.=FALSE)
} else {
  checkp(ks.test(Fx1(X1$x), "punif")$p.value,
       the.context,
       "Kolmogorov-Smirnov test of uniformity of transformed x coordinates of type 1 points")
}
checkp(ks.test(X2$x, "punif")$p.value,
       the.context,
     "Kolmogorov-Smirnov test of uniformity of x coordinates of type 2 points")
checkp(ks.test(Fy2(X2$y), "punif")$p.value,
       the.context,
       "Kolmogorov-Smirnov test of uniformity of transformed y coordinates of type 2 points")
