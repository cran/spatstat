# test a scoping problem that can arise when ppm splits the data

require(spatstat)
data(bei)
fit <- ppm(bei, ~elev, covariates=bei.extra)
mm <- model.matrix(fit)
