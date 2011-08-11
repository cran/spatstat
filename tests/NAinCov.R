# tests/NAinCov.R
# Testing the response to the presence of NA's in covariates

require(spatstat)
X <- runifpoint(42)
Y <- as.im(function(x,y) { x+y }, owin())
Y[owin(c(0.2,0.4),c(0.2,0.4))] <- NA
# fit model: should produce a warning but no failure
misfit <- ppm(X, ~Y, covariates=list(Y=Y))
# prediction 
Z <- predict(misfit, type="trend")
Z <- predict(misfit, type="se")
# covariance matrix: all should be silent
v <- vcov(misfit)
ss <- vcov(misfit, what="internals")





