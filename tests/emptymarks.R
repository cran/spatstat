#
# test cases where there are no (rows or columns of) marks
#

require(spatstat)
data(cells)
n <- npoints(cells)
df <- data.frame(x=1:n, y=factor(sample(letters, n, replace=TRUE)))
nocolumns <- c(FALSE, FALSE)
norows <- rep(FALSE, n)
X <- cells
marks(X) <- df
marks(X) <- df[,1]
marks(X) <- df[,nocolumns]
Z <- Y <- X[integer(0)]
marks(Y) <- df[norows,]
stopifnot(is.marked(Y))
marks(Z) <- df[norows,nocolumns]
stopifnot(!is.marked(Z))


