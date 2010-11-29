clarkevans <- function(X, correction=c("none", "Donnelly", "cdf"),
                       clipregion=NULL)
{
  verifyclass(X, "ppp")
  W <- X$window
  area <- area.owin(W)
  npts <- npoints(X)
  intensity <- npts/area

  # R undefined for empty point pattern
  if(npts == 0)
    return(NA)

  # validate correction argument
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Donnelly="Donnelly",
                             donnelly="Donnelly",
                             guard="guard",
                             cdf="cdf"),
                           multi=TRUE)

  if(("Donnelly" %in% correction) && (W$type != "rectangle"))
    warning("Donnelly correction only available for rectangular windows")

  # guard correction applied iff `clipregion' is present
  askguard <- ("guard" %in% correction)
  gaveguard <- !is.null(clipregion)
  if(askguard && !gaveguard)
    warning("guard correction not performed; clipregion not specified")
  else if(gaveguard && !askguard)
    correction <- c(correction, "guard")
  
  # Dobs = observed mean nearest neighbour distance
  nndistX <- nndist(X)
  Dobs <- mean(nndistX)
  # Dpois = Expected mean nearest neighbour distance for Poisson process
  Dpois <- 1/(2*sqrt(intensity))

  answer <- NULL
  
  # Naive uncorrected value
  if("none" %in% correction) {
    Rnaive <- Dobs/Dpois
    answer <- c(answer, naive=Rnaive)
  }
  # Donnelly edge correction
  if("Donnelly" %in% correction) {
     # Dedge = Edge corrected mean nearest neighbour distance, Donnelly 1978
    if(W$type == "rectangle") {
      perimeter <- 2*(diff(W$xrange) + diff(W$yrange))
      Dkevin  <- Dpois + (0.0514+0.0412/sqrt(npts))*perimeter/npts
      Rkevin <- Dobs/Dkevin
    } else 
      Rkevin <- NA
    answer <- c(answer, Donnelly=Rkevin)
  }
  # guard area method
  if("guard" %in% correction && !is.null(clipregion)) {
    # use nn distances from points inside `clipregion'
    clip <- as.owin(clipregion)
    ok <- inside.owin(X, , clip)
    Dguard <- mean(nndistX[ok])
    Rguard <- Dguard/Dpois
    answer <- c(answer, guard=Rguard)
  }
  if("cdf" %in% correction) {
    # compute mean of estimated nearest-neighbour distance distribution G
    G <- Gest(X)
    numer <- stieltjes(function(x){x}, G)$km
    denom <- stieltjes(function(x){rep(1, length(x))}, G)$km
    Dcdf <- numer/denom
    Rcdf <- Dcdf/Dpois
    answer <- c(answer, cdf=Rcdf)
  }
  return(answer)
}

