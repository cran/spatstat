require(spatstat, quietly=TRUE, save=FALSE)
.getLetterR <- function() {
  filename <- paste(spatstat.rawdata.location(), "letterR.tab",
                    sep=.Platform$file.sep)
  lt <- read.table(filename, header=TRUE)
  hole <- lt[lt$which == "inner", ]
  outside <- lt[lt$which == "outer", ]
  letterR <- owin(poly=list(
                  list(x=outside$x,y=outside$y),
                  list(x=hole$x,y=hole$y)))
  letterR <- affine(letterR, mat=diag(c(1,1)/1000))
  return(letterR)
}
letterR <- .getLetterR()
rm(.getLetterR)


