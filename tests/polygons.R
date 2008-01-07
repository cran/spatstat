require(spatstat)
data(letterR)
co <- as.ppp(corners(letterR), letterR, check=FALSE)
co[letterR]
