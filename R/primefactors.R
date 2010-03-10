#
#  primefactors.R
#
#  $Revision: 1.3 $   $Date: 2005/05/10 01:32:36 $
#

eratosthenes <- function(nmax) {
  if(nmax < 2) return(numeric(0))
  numbers <- 2:nmax
  prime <- 2
  repeat{
    retain <-  (numbers <= prime) | (numbers %% prime != 0)
    numbers <- numbers[retain]
    remaining <- (numbers > prime)
    if(!any(remaining))
      break
    prime <- min(numbers[remaining])
  }
  return(numbers)
}
  
primefactors <- function(n, prmax) {
  if(missing(prmax)) prmax <- floor(sqrt(n))
  primes <- eratosthenes(prmax)
  divides.n <- (n %% primes == 0)
  if(!any(divides.n)) 
    return(n)
  else {
    divisors <- primes[divides.n]
    prmax <- max(divisors)
    m <- n/prod(divisors)
    if(m == 1) return(divisors)
    else {
      mfactors <- primefactors(m, prmax=prmax)
      return(sort(c(divisors, mfactors)))
    }
  }
}

is.prime <- function(n) { length(primefactors(n)) == 1 }

least.common.multiple <- function(n, m) {
  nf <- primefactors(n)
  mf <- primefactors(m)
  p <- sort(unique(c(nf,mf)))
  nfac <- table(factor(nf, levels=p))
  mfac <- table(factor(mf, levels=p))
  prod(p^pmax(nfac,mfac))
}

greatest.common.divisor <- function(n, m) {
  nf <- primefactors(n)
  mf <- primefactors(m)
  p <- sort(unique(c(nf,mf)))
  nfac <- table(factor(nf, levels=p))
  mfac <- table(factor(mf, levels=p))
  prod(p^pmin(nfac,mfac))
}
  
divisors <- function(n) {
  p <- primefactors(n)

  up <- sort(unique(p))
  k <- table(factor(p, levels=up))

  rest <- function(kk, uu) {
    powers <- uu[1]^(0:(kk[1]))
    if(length(uu) == 1)
      return(powers)
    rr <- rest(kk[-1], uu[-1])
    products <- as.vector(outer(powers, rr, "*"))
    return(sort(unique(products)))
    }

  return(rest(k, up))
}

    
