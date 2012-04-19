#
# Asymptotic covariance & correlation matrices
# and Fisher information matrix
# for ppm objects
#
#  $Revision: 1.32 $  $Date: 2012/04/15 01:03:32 $
#

vcov.ppm <- local({

vcov.ppm <- function(object, ..., what="vcov", verbose=TRUE,
                     gam.action=c("warn", "fatal", "silent"),
                     matrix.action=c("warn", "fatal", "silent"),
                     hessian=FALSE) {
  verifyclass(object, "ppm")

  gam.action <- match.arg(gam.action)
  matrix.action <- match.arg(matrix.action)

  stopifnot(length(what) == 1 && is.character(what))
  what.options <- c("vcov", "corr", "fisher", "Fisher", "internals")
  what.map     <- c("vcov", "corr", "fisher", "fisher", "internals")
  if(is.na(m <- pmatch(what, what.options)))
    stop(paste("Unrecognised option: what=", sQuote(what)))
  what <- what.map[m]

  method <- resolve.defaults(list(...), list(method="C"))$method
  
  # Fisher information *may* be contained in object
  fisher <- object$fisher
  varcov <- object$varcov
  
  # Do we need to go into the guts?
  needguts <-
    (is.null(fisher) && what=="fisher") ||
    (is.null(varcov) && what %in% c("vcov", "corr")) ||
    (what == "internals")

  # In general it is not true that varcov = solve(fisher)
  # because we might use different estimators,
  # or the parameters might be a subset of the canonical parameter

  ############## guts ##############################

  if(needguts) {

    # warn if fitted model was obtained using GAM
    if(identical(object$fitter, "gam")) {
      switch(gam.action,
             fatal={
               stop(paste("model was fitted by gam();",
                          "execution halted because fatal=TRUE"),
                    call.=FALSE)
             },
             warn={
               warning(paste("model was fitted by gam();",
                             "asymptotic variance calculation ignores this"),
                       call.=FALSE)
             },
             silent={})
    }

    if(!is.poisson(object) && !hessian) {
      # Gibbs model - invoke new code from Ege and J-F
      varcov <- vcovGibbs(object, ..., matrix.action=matrix.action)
      if(is.null(varcov))
        return(NULL)
      if(what %in% c("fisher", "internals")) {
        # compute fisher information 
        fisher <- checksolve(varcov, matrix.action,
                             "variance-covariance matrix",
                             "Fisher information" )
        if(is.null(fisher)) return(NULL)
      }
    } else {
      # Poisson model, or Hessian of Gibbs model
      gf <- getglmfit(object)
      # we need a glm or gam
      if(is.null(gf)) {
        if(verbose) 
          warning("Refitting the model using GLM/GAM")
        object <- update(object, forcefit=TRUE)
        gf <- getglmfit(object)
        if(is.null(gf))
          stop("Internal error - refitting did not yield a glm object")
      }
      # compute related stuff
      fi <- fitted(object, type="trend", check=FALSE, drop=TRUE)
      wt <- quad.ppm(object, drop=TRUE)$w
      # extract model matrix
      mom <- model.matrix(object, keepNA=FALSE) 
      # compute Fisher information if not known
      if(is.null(fisher)) {
        switch(method,
               C = {
                 fisher <- sumouter(mom, fi * wt)
               },
               interpreted = {
                 for(i in 1:nrow(mom)) {
                   ro <- mom[i, ]
                   v <- outer(ro, ro, "*") * fi[i] * wt[i]
                   if(!any(is.na(v)))
                     fisher <- fisher + v
                 }
                 momnames <- dimnames(mom)[[2]]
                 dimnames(fisher) <- list(momnames, momnames)
               })
      }
      
    }

  }
  
  ############## end of guts ####################################
  
  if(what %in% c("vcov", "corr") && is.null(varcov)) {
    # Need variance-covariance matrix.
    # Derive from Fisher information
    varcov <- checksolve(fisher, matrix.action,
                         "Fisher information matrix",
                         "variance")
    if(is.null(varcov)) return(NULL)
  }
         
  switch(what,
         fisher = { return(fisher) },
         vcov   = {
           return(varcov)
         },
         corr={
           sd <- sqrt(diag(varcov))
           return(varcov / outer(sd, sd, "*"))
         },
         internals = {
           return(list(fisher=fisher, suff=mom))
         })
}

  # vcovGibbs from Ege Rubak and J-F Coeurjolly

vcovGibbs <- function(fit, ..., special.alg = TRUE,
                      matrix.action=c("warn", "fatal", "silent")){
  verifyclass(fit, "ppm")
  if(!is.stationary(fit))
    stop("Sorry, vcov is only implemented for stationary processes")
  if(fit$correction!="border")
    stop("Sorry, vcov is only implemented for the border correction")

  matrix.action <- match.arg(matrix.action)
  
  ## Interaction name:
  iname <- fit$interaction$name
  
  ## Handling the Area Interaction Process: no longer needed
#  if(iname=="Area-interaction process")
#    fit$interaction$pot <- HackedAreaInterPot
  
  ## Does the model have marks which are in the trend?
  marx <- is.marked(fit) && ("marks" %in% variablesinformula(fit$trend))

  ## The full data and window:
  Xplus <- data.ppm(fit)
  Wplus <- as.owin(Xplus)

  ## Fitted parameters and the parameter dimension p (later consiting of p1 trend param. and p2 interaction param.):
  theta <- coef(fit)
  p <- length(theta)
  
  ## Number of points:
  n <- npoints(Xplus)

  ## Using the faster algorithms for special cases
  if(special.alg){
    param <- coef(fit)
    switch(iname,
      "Strauss"={
	  return(vcovPairPiece(Xplus, reach(fit$interaction), exp(coef(fit)[2])))
      },
           
      "Piecewise constant pairwise interaction process"={
        return(vcovPairPiece(Xplus, fit$interaction$par$r, exp(coef(fit)[-1])))
      },

      "Multitype Strauss process"={
	matR <- fit$interaction$par$radii
        R <- c(matR[1,1], matR[1,2], matR[2,2])
        ## Only implemented for 2 types with equal interaction range:
        if(ncol(matR)==2 && marx){
          n <- length(theta)
          res <- vcovMultiStrauss(Xplus, R, exp(theta[c(n-2,n-1,n)]))
          return(res)
        }
      }
    )
  }
  
  ## Matrix specifying equal points in the two patterns in the call to eval below:
  E <- matrix(rep(1:n, 2), ncol = 2)

  ## Eval. the interaction potential difference at all points (internal spatstat function):
  V1 <- fit$interaction$family$eval(Xplus, Xplus, E, fit$interaction$pot, fit$interaction$par, fit$correction)

  ## Calculate parameter dimensions and correct the contrast type parameters:
  p2 <- ncol(V1)
  p1 <- p-p2
  if(p1>1)
    theta[2:p1] <- theta[2:p1] + theta[1]
  ## V1 <- evalInteraction(Q, Xplus, union.quad(Q), fit$interaction, fit$correction)
  POT <- attr(V1, "POT")
  attr(V1, "POT") <- NULL
  ## Adding the constant potential as first column (one column per type for multitype):
  if(!marx){
    V1 <- cbind(1, V1)
    colnames(V1) <- names(theta)
  }
  else{
    lev <- levels(marks(Xplus))
    ## Indicator matrix for mark type attached to V1:
    tmp <- matrix(marks(Xplus), nrow(V1), p1)==matrix(lev, nrow(V1), p-ncol(V1), byrow=TRUE)
    colnames(tmp) <- lev
    V1 <- cbind(tmp,V1)
  }

  ## Matrices for differences of potentials:
  E <- matrix(rep(1:(n-1), 2), ncol = 2)
  V2 <- array(0,dim=c(n,n,p))
  dV <- array(0,dim=c(n,n,p))

  for(k in 1:p1){
    V2[,,k] <- matrix(V1[,k], n, n, byrow = FALSE)
  }
  for(k in (p1+1):p){
    diag(V2[,,k]) <- V1[,k]
  }
  for(j in 1:n){
    ## Fast evaluation for pairwise interaction processes:
    if(fit$interaction$family$name=="pairwise"){
      V2[-j,j,-(1:p1)] <- V1[-j,-(1:p1)]-POT[-j,j,]
    }
    else{
      V2[-j,j,-(1:p1)] <- fit$interaction$family$eval(Xplus[-j], Xplus[-j], E, fit$interaction$pot, fit$interaction$par, fit$correction)
      ## Q <- quadscheme(Xplus[-j],emptyppp)
      ## V2[-j,j,-1] <- evalInteraction(Q, Xplus[-j], Xplus[-j], fit$interaction, fit$correction)
    }
    for(k in 1:p){
      dV[,j,k] <- V1[,k] - V2[,j,k]
    }
  }
  ## Ratio of first and second order Papangelou - 1:
  frac <- 0*dV[,,1]
  for(k in (p1+1):p){
    frac <- frac + dV[,,k]*theta[k]
  }
  frac <- exp(-frac)-1

  ## In the rest we restric attention to points in the interior:
  
  ## The interaction range:
  R <- reach(fit$interaction)

  ## The reduced window, area and point pattern:
  W<-erosion.owin(Wplus,R)
  areaW <- area.owin(W)
  X <- Xplus[W]
  
  ## Making a logical matrix, I, indicating R-close pairs which are in the interior:
  IntPoints <- nncross(Xplus,X)$dist==0
  D <- pairdist(Xplus)
  diag(D) <- Inf
  I <- D<=R * outer(IntPoints,IntPoints)
  
  ## Matrix A1:
  A1 <- t(V1[IntPoints,])%*%V1[IntPoints,]/areaW

  ## Matrix A2:
  A2 <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A2[k,l] <- A2[l,k] <- sum(I*V2[,,k]*frac*t(V2[,,l]))
    }
  }
  A2 <- A2/areaW
  
  ## Matrix A3:
  A3 <- matrix(0,p,p)
  for(k in 1:p){
    for(l in k:p){
      A3[k,l] <- A3[l,k] <- sum(I*dV[,,k]*t(dV[,,l]))
    }
  }
  A3 <- A3/areaW
  
  ## Matrix Sigma (A1+A2+A3):
  Sigma<-A1+A2+A3
  
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(NULL)

  mat <- U%*%Sigma%*%U / areaW

  ## Convert to treatment contrasts
 if(marx){
    mat <- contrastmatrix(mat, p1)
    dimnames(mat) <- list(names(theta), names(theta))
  }
  
  return(mat)
}

vcovPairPiece <- function(Xplus, R, gam, matrix.action){
  ## R is  the  vector of breaks (R[length(R)]= range of the pp.
  ## gam is the vector of weights
  
  ## Xplus : point process observed in W+R
  ## Extracting the window and calculating area:
  Wplus<-as.owin(Xplus)
  W<-erosion.owin(Wplus,R[length(R)])
  areaW <- area.owin(W)
  X<-Xplus[W]
  # Identify points inside W
  IntPoints<- inside.owin(Xplus, , W)
  # IntPoints <- (Xplus$x >= W$xrange[1]) &  (Xplus$x <= W$xrange[2])&(Xplus$y >= W$yrange[1]) &  (Xplus$y <= W$yrange[2])
  
  ## Matrix D with pairwise distances between points and infinite distance
  ## between a point and itself:
  
  Dplus<-pairdist(Xplus)
  D <- pairdist(X)
  diag(D) <- diag(Dplus) <- Inf
  ## logical matrix, I, indicating R-close pairs:
  p<-length(R)
  Tplus<-T<-matrix(0,X$n,p)
  I<-Iplus<-list()
  for (i in 1:p){
     if (i==1){
	Iplus[[1]]<- Dplus <=R[1]
	I[[1]] <- D<=R[1]
     } else {
	Iplus[[i]]<- ((Dplus>R[i-1]) & (Dplus <=R[i]))
	I[[i]] <- ((D>R[i-1]) & (D <=R[i]))
     }
     ## Vector T with the number of $R$-close neighbours to each point:
     Tplus[,i]<-colSums(Iplus[[i]])[IntPoints]
     T[,i] <- colSums(I[[i]])
  }
  ## Matrices A1, A2 and A3 are initialized to zero:
  A1 <- A2 <- A3 <- matrix(0,p+1,p+1)
  ## A1 and A3:
  A1[1,1] <- npoints(X)
  
  for (j in (2:(p+1))){
    A1[1,j]<-A1[j,1]<-sum(Tplus[,j-1])
    A3[j,j]<-sum(T[,j-1])
    for (k in (2:(p+1))){
      A1[j,k]<-sum(Tplus[,j-1] * Tplus[,k-1])
    }
  }
  ## A2:
  for (j in (2:(p+1))){
    A2[1,1]<-A2[1,1]+(gam[j-1]^(-1)-1)*sum(T[,j-1])
    for (l in (2:(p+1))){
      if (l==j) vj<-Tplus[,j-1]-1 else vj<-Tplus[,j-1]
	A2[1,j]<-A2[1,j]+(gam[l-1]^(-1)-1)*sum(T[,l-1]*(vj) )
    }
    A2[j,1]<-A2[1,j]
    for (k in (2:(p+1))){
      for (l in (2:(p+1))){
	if (l==j) vj<-Tplus[,j-1]-1 else vj<-Tplus[,j-1]
	if (l==k) vk<-Tplus[,k-1]-1 else vk<-Tplus[,k-1]

	A2[j,k]<-A2[j,k]+ (gam[l-1]^(-1)-1)*sum(I[[l-1]]*outer(vj,vk))
      }
    }

  }
  A1<-A1/areaW
  A2<-A2/areaW
  A3<-A3/areaW
  ## browser()
  Sigma<-A1+A2+A3
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(NULL)

  mat <- U%*%Sigma%*%U / areaW
  nam <- c("(intercept)", names(gam))
  dimnames(mat) <- list(nam, nam)
  return(mat)
}

vcovMultiStrauss <- function(Xplus, vecR, vecg, matrix.action){
  ## Xplus : marked Strauss point process 
  ## with two types 
  ## observed in W+R (R=max(R11,R12,R22))

  ## vecg =  estimated parameters of interaction parameters
  ##	    ordered as the output of ppm, i.e. vecg=(g11,g12,g22)	
  ## vecR = range for the diff. strauss ordered a vecg(R11,R12,R22)

  R <- max(vecR)
  R11<-vecR[1];R12<-vecR[2];R22<-vecR[3]
  ## Extracting the window and calculating area:
  Wplus<-as.owin(Xplus)
  W<-erosion.owin(Wplus,R)
  areaW <- area.owin(W)
  X<-Xplus[W]
  X1plus<-Xplus[Xplus$marks==levels(Xplus$marks)[1]]
  X2plus<-Xplus[Xplus$marks==levels(Xplus$marks)[2]]
  X1<-X1plus[W]
  X2<-X2plus[W]
  # Identify points inside W
  IntPoints1 <- inside.owin(X1plus, , W)
  IntPoints2 <- inside.owin(X2plus, , W)
  # IntPoints1<- (X1plus$x >= W$xrange[1]) &  (X1plus$x <= W$xrange[2])&(X1plus$y >= W$yrange[1]) &  (X1plus$y <= W$yrange[2]) 
  # IntPoints2<- (X2plus$x >= W$xrange[1]) &  (X2plus$x <= W$xrange[2])&(X2plus$y >= W$yrange[1]) &  (X2plus$y <= W$yrange[2]) 

  ## Matrix D with pairwise distances between points and infinite distance
  ## between a point and itself:

  D1plus<-pairdist(X1plus)
  D1 <- pairdist(X1)
  diag(D1) <- diag(D1plus) <- Inf
  
  D2plus<-pairdist(X2plus)
  D2 <- pairdist(X2)
  diag(D2) <- diag(D2plus) <- Inf
  
  D12plus<-crossdist(X1,X2plus)  
  T12plus<-rowSums(D12plus<=R12)
  D21plus<-crossdist(X2,X1plus) 
  T21plus<-rowSums(D21plus<=R12)
  
  I12<-crossdist(X1,X2)<=R12
  I21<-crossdist(X2,X1)<=R12
  T12<-rowSums( I12)  
  T21<-rowSums(I21)
  ## logical matrix, I, indicating R-close pairs:
  I1plus<- D1plus <=R11
  I1 <- D1<=R11
  I2plus<- D2plus <=R22
  I2 <- D2<=R22
  ## Vector T with the number of $R$-close neighbours to each point:
  T1plus<-colSums(I1plus)[IntPoints1]
  T1 <- colSums(I1)
  T2plus<-colSums(I2plus)[IntPoints2]
  T2 <- colSums(I2)

  ## Matrices A1, A2 and A3 are initialized to zero:
  A1 <- A2 <- A3 <- matrix(0,5,5)
  ## A1 is filled:
  A1[1,1]<-npoints(X1)
  A1[1,3]<-A1[3,1]<-sum(T1plus)
  A1[1,4]<-A1[4,1]<-sum(T12plus)
  A1[2,2]<-npoints(X2)
  A1[2,5]<-A1[5,2]<-sum(T2plus)
  A1[2,4]<-A1[4,2]<-sum(T21plus)
  A1[3,3]<-sum(T1plus*T1plus)
  A1[3,4]<-A1[4,3]<-sum(T1plus*T12plus)
  A1[5,5]<-sum(T2plus*T2plus)
  A1[4,5]<-A1[5,4]<-sum(T2plus*T21plus)
  A1[4,4]<-sum(T12plus*T12plus)+sum(T21plus*T21plus)

  ## A3 is filled:
  A3[3,3]<-sum(T1)
  A3[5,5]<-sum(T2)
  A3[4,4]<-sum(T12)+sum(T21)
   

  ## A2 is filled:
  gamInv<-vecg^(-1)-1
  gi1<-gamInv[1];gi12<-gamInv[2];gi2<-gamInv[3]
  A2[1,1]<-sum(T1)*gi1
  A2[1,2]<-A2[2,1]<-sum(T12)*gi12
  A2[1,3]<-A2[3,1]<-sum(T1*(T1plus-1))*gi1
  A2[1,5]<-A2[5,1]<-sum(T21*T2plus)*gi12
  A2[1,4]<-A2[4,1]<-gi1*sum(T1*(T12plus))+gi12*sum(T21*(T21plus-1))
  A2[2,2]<-sum(T2)*gi2
  A2[2,3]<-A2[3,2]<-sum(T12*T1plus)*gi12
  A2[2,5]<-A2[5,2]<-sum(T2*(T2plus-1))*gi2
  A2[2,4]<-A2[4,2]<-gi2*sum(T2*(T21plus))+gi12*sum(T12*(T12plus-1))

  A2[3,3]<-gi1*sum(I1*outer(T1plus-1,T1plus-1))
  
  A2[3,5]<-A2[5,3]<- gi12*sum(I12*outer(T1plus,T2plus))
  A2[3,4]<-A2[4,3]<-gi1*sum(I1*outer(T1plus-1,T12plus))+gi12*sum(I12*outer(T1plus,T21plus-1))
  
  A2[5,5]<-gi2*sum(I2*outer(T2plus-1,T2plus-1))
  A2[4,5]<-A2[5,4]<-gi2*sum(I2*outer(T2plus-1,T21plus))+gi12*sum(I21*outer(T2plus,T12plus-1))
  
  A2[4,4]<-gi1*sum(I1*outer(T12plus,T12plus))+gi2*sum(I2*outer(T21plus,T21plus))+ gi12*sum(I12*outer(T12plus-1,T21plus-1))+gi12*sum(I21*outer(T21plus-1,T12plus-1))
  
  #browser()
  A1<-A1/areaW
  A2<-A2/areaW
  A3<-A3/areaW
  #browser()
  
  Sigma<-A1+A2+A3
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(NULL)

  mat <- U%*%Sigma%*%U / areaW

  nam <- c(levels(marks(Xplus)), names(vecg))
  dimnames(mat) <- list(nam, nam)
  
  return(mat)
}

#############################################################
#############  only for strauss...included in vcovPairPiece!
#############################################################


vcovStrauss <- function(Xplus, R, gam, matrix.action){
  ## Xplus : point process observed in W+R
  ## Extracting the window and calculating area:
  Wplus<-as.owin(Xplus)
  W<-erosion.owin(Wplus,R)
  areaW <- area.owin(W)
  X<-Xplus[W]
  # Identify points inside W
  # IntPoints<- (Xplus$x >= W$xrange[1]) &  (Xplus$x <= W$xrange[2])&(Xplus$y >= W$yrange[1]) &  (Xplus$y <= W$yrange[2])
  IntPoints <- inside.owin(Xplus, , W)
  
  ## Matrix D with pairwise distances between points and infinite distance
  ## between a point and itself:
  
  Dplus<-pairdist(Xplus)
  D <- pairdist(X)
  diag(D) <- diag(Dplus) <- Inf
  ## logical matrix, I, indicating R-close pairs:
  Iplus<- Dplus <=R
  I <- D<=R
  ## Vector T with the number of $R$-close neighbours to each point:
  Tplus<-colSums(Iplus)[IntPoints]
  T <- colSums(I)
  
  ## Matrices A1, A2 and A3 are initialized to zero:
  A1 <- A2 <- A3 <- matrix(0,2,2)
  ## A1 is filled:
  A1[1,1] <- npoints(X)
  A1[1,2] <- A1[2,1] <- sum(Tplus)
  A1[2,2] <- sum(Tplus*Tplus)
  A1 <- A1/areaW
  if (A1[2,2]==0) A1[2,2]<-A1[1,2]<-A1[2,1]<-1
  ## A2 is filled:
  A2[1,1] <- sum(T)
  A2[1,2] <- A2[2,1] <- sum(T*(Tplus-1))
  A2[2,2] <- sum(I*outer(Tplus-1,Tplus-1))
  A2 <- A2*(gam^(-1)-1)/areaW
  ## A3 is filled:
  A3[2,2] <- sum(T)/areaW
  Sigma<-A1+A2+A3
  ## Finally the result is calculated:
  U <- checksolve(A1, matrix.action, , "variance")
  if(is.null(U)) return(U)

  mat <- U%*%Sigma%*%U / areaW
  
  return(mat)
}

# Convert the first p rows & columns of variance matrix x
# to variances of treatment contrasts
contrastmatrix <- function(x,p){
  mat <- x
  ## Correct column and row 1:
  for(i in 2:p){
    mat[1,i] <- mat[i,1] <- x[1,1]-x[1,i]
  }
  ## Correct columns and rows 2,...,p:
  for(i in 2:p){
    for(j in 2:p){
      mat[i,j] <- x[1,1]-x[1,i]-x[1,j]+x[i,j]
    }
    for(j in (p+1):ncol(x)){
      mat[i,j] <- mat[j,i] <- x[1,j]-x[i,j]
    }
  }
  mat
}


checksolve <- function(M, action, descrip, target="") {
  Mname <- deparse(substitute(M))
  Minv <- try(solve(M), silent=(action=="silent"))
  if(!inherits(Minv, "try-error"))
    return(Minv)
  if(missing(descrip))
    descrip <- paste("the matrix", Mname)
  whinge <- paste("Cannot compute",
                  paste(target, ":", sep=""),
                  descrip, "is singular")
  switch(action,
         fatal=stop(whinge, call.=FALSE),
         warn= warning(whinge, call.=FALSE),
         silent={})
  return(NULL)
}

vcov.ppm
}
)

suffloc <- function(object) {
  verifyclass(object, "ppm")
  if(!is.poisson(object))
    stop("Internals not available for Gibbs models")
  return(vcov(object, what="internals")$suff)
}

