#
# Determine which 'canonical variables' depend on a supplied covariate
#
#   $Revision: 1.2 $  $Date: 2010/05/25 11:04:45 $
#

model.depends <- function(object) {
  # supplied covariates
  fo <- formula(object)
  if(length(as.list(fo)) == 3) {
    # formula has a response: strip it
    fo <- fo[-2]
  }
  covars <- variablesinformula(fo)
  # canonical covariates 
  mm <- model.matrix(object)
  ass <- attr(mm, "assign")
  # model terms
  tt <- terms(object)
  lab <- attr(tt, "term.labels")
  # 'ass' maps canonical covariates to 'lab'
  # determine which canonical covariate depends on which supplied covariate
  depends <- matrix(FALSE, length(ass), length(covars))
  for(i in seq(along=ass)) {
    if(ass[i] == 0) # 0 is the intercept term
      depends[i,] <- FALSE
    else {
      turm <- lab[ass[i]]
      depends[i, ] <- covars %in% all.vars(parse(text=turm))
    }
  }
  rownames(depends) <- colnames(mm)
  colnames(depends) <- covars
  return(depends)
}

model.is.additive <- function(object) {
  dep <- model.depends(object)
  hit <- t(dep) %*% dep
  diag(hit) <- 0
  ok <- all(hit == 0)
  return(ok)
}
