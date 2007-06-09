clarkevans <- function(X, correction=c("none", "Donnelly", "guard"),
                       clipregion=NULL)
{
  verifyclass(X, "ppp")
  W <- X$window
  area <- area.owin(W)
  npoints <- X$n
  intensity <- npoints/area

  # R undefined for empty point pattern
  if(npoints == 0)
    return(NA)

  # validate correction argument
  correction <- pickoption("correction", correction,
                           c(none="none",
                             Donnelly="Donnelly",
                             donnelly="Donnelly",
                             guard="guard"),
                           multi=TRUE)

  if(("Donnelly" %in% correction) && (W$type != "rectangle"))
    warning("Donnelly correction only available for rectangular windows")
    
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
      Dedge  <- Dpois + (0.0514+0.0412/sqrt(npoints))*perimeter/npoints
      Redge <- Dobs/Dedge
    } else 
      Redge <- NA
    answer <- c(answer, edge=Redge)
  }
  # guard area method
  if("guard" %in% correction) {
    if(!is.null(clipregion)) {
      # use nn distances from points inside `clipregion'
      clip <- as.owin(clipregion)
      ok <- inside.owin(X$x, X$y, clip)
    } else {
      # use uncensored nn distances
      bdistX <- bdist.points(X)
      ok <- (nndistX < bdistX)
    }
    Dguard <- mean(nndistX[ok])
    Rguard <- Dguard/Dpois
    answer <- c(answer, guard=Rguard)
  }
  return(answer)
}

