# Berman-Huntington points & lines data
                             
.spatstat.get.copper <- function() {
  copper <- list()
  copper$FullWindow <- owin(c(-0.335, 70.11), c(0.19, 158.233))

  # read data
  pointsfile <- spatstat.rawdata.location("copper-points.tab")
  linesfile <- spatstat.rawdata.location("copper-lines.tab")
  df <- read.table(pointsfile, header=TRUE)
  copper$points <- ppp(df$x, df$y, window=copper$FullWindow)
  copper$lines <- read.table(linesfile, header=TRUE)[,2:5]

  # clip to Southern half window
  copper$SouthWindow <- owin(c(-0.335, 35), c(0.19, 158.233))
  copper$SouthPoints <- copper$points[, copper$SouthWindow]
  copper$SouthLines <- .spatstat.clip35(copper$lines)

  # code to compute distance to nearest lineament
  copper$SouthDistance <- function() {
    bw <- as.mask(copper$SouthWindow)
    xx <- raster.x(bw)
    yy <- raster.y(bw)
    dd <- distppll(cbind(as.vector(xx),as.vector(yy)), copper$lines)
    dmin <- apply(dd, 1, min)
    result <- im(matrix(dmin, ncol=ncol(bw$m), nrow=nrow(bw$m)),
                      xcol=bw$xcol, yrow=bw$yrow)
    return(result)
  }

  return(copper)
}

.spatstat.clip35 <- function(segs) {
# original author: Rob Foxall 1997
# hacked by Adrian Baddeley 2004

# The function clip35() clips the line segments
# to the subregion where x <= 35.
        xl1 <- segs$x1
        yl1 <- segs$y1
        xl2 <- segs$x2
        yl2 <- segs$y2
        
	nlines <- length(xl1)
	xmax <- rep(0, nlines)
	xmin <- xmax
	for(i in 1:nlines) {
		xmax[i] <- max(xl1[i], xl2[i])
		xmin[i] <- min(xl1[i], xl2[i])
	}
	xlok <- xmin <= 35
	n <- sum(xlok)
	x1 <- rep(0, n)
	y1 <- x1
	x2 <- x1
	y2 <- x1
	i <- 1
	for(j in 1:nlines) {
		if(xmax[j] <= 35) {
			x1[i] <- xl1[j]
			y1[i] <- yl1[j]
			x2[i] <- xl2[j]
			y2[i] <- yl2[j]
			i <- i + 1
		}
		if(xmax[j] >= 35 & xmin[j] <= 35) {
			if(xl1[j] <= 35) {
				x1[i] <- xl1[j]
				y1[i] <- yl1[j]
				x2[i] <- 35
				y2[i] <- yl1[j] + ((35 - xl1[j])/(xl2[j] - 
				  xl1[j])) * (yl2[j] - yl1[j])
				i <- i + 1
			}
			if(xl2[j] <= 35) {
				x1[i] <- xl2[j]
				y1[i] <- yl2[j]
				x2[i] <- 35
				y2[i] <- yl2[j] + ((35 - xl2[j])/(xl1[j] - 
				  xl2[j])) * (yl1[j] - yl2[j])
				i <- i + 1
			}
		}
	}
	return(data.frame(x1=x1, y1=y1, x2=x2, y2=y2))
}

########### do it ###################

copper <- .spatstat.get.copper()

# clean up
rm(.spatstat.get.copper, .spatstat.clip35)


