# test of case where mark levels contain illegal characters

require(spatstat)
hyphenated <- c("a", "not-a")
spaced <- c("U", "non U")
suffixed <- c("a+", "a*")
charred <- c("+", "*")

data(amacrine)
irad <- matrix(0.1, 2,2)
hrad <- matrix(0.005, 2, 2)

tryit <- function(types, X, irad, hrad) { 
  m <- X$marks
  levels(m) <- types
  X$marks <- m
  ppm(X, ~marks + polynom(x,y,2), MultiStrauss(types=types,radii=irad))
  fit <- ppm(X, ~marks + polynom(x,y,2),
    MultiStraussHard(types=types,iradii=irad,hradii=hrad))
  print(fit)
  print(coef(fit))
  val <- fitted(fit)
  pred <- predict(fit)
  return(NULL)
}


tryit(hyphenated, amacrine, irad, hrad)
tryit(spaced, amacrine, irad, hrad)
tryit(suffixed, amacrine, irad, hrad)
tryit(charred, amacrine, irad, hrad)
