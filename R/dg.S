#
#     dg.S
#
#     Experimental
#
#     Diggle-Gratton pair potential
#
#
dgpot <- function(tt, par) {
  delta <- par$delta
  rho <- par$rho
  above <- (tt > rho)
  inrange <- (!above) & (tt > delta)
  h <- above + inrange * (tt - delta)/(rho - delta)
  log(h)
}

DiggleGratton <- function(delta, rho) {
  Pairwise(dgpot, "Diggle-Gratton process",
           list(delta=delta,rho=rho),
           list("lower limit delta", "upper limit rho"), NULL)
}

