#
# hypersub.R
#
#
#  subset operations for hyperframes
#
#  $Revision: 1.4 $    $Date: 2010/01/14 19:40:17 $
#
#

"[.hyperframe" <- function(x, i, j, drop=FALSE, ...) {
  x <- unclass(x)
  if(!missing(i)) {
    y <- x
    y$df     <- x$df[i, , drop=FALSE]
    y$ncases <- nrow(y$df)
    y$hypercolumns <- lapply(x$hypercolumns, function(z,k) { z[k] }, k=i)
    x <- y
  }
  if(!missing(j)) {
    y <- x
    patsy <- seq(y$nvars)
    names(patsy) <- y$vname
    jj <- patsy[j]
    names(jj) <- NULL
    y$nvars <- length(jj)
    y$vname <- vname <- x$vname[jj]
    y$vtype <- vtype <- x$vtype[jj]
    y$vclass <- x$vclass[jj]
    if(ncol(x$df) != 0) 
      y$df    <- x$df[ , vname[vtype == "dfcolumn"], drop=FALSE]
    y$hyperatoms <- x$hyperatoms[ vname[ vtype == "hyperatom" ]]
    y$hypercolumns <- x$hypercolumns[ vname [ vtype == "hypercolumn" ] ]
    x <- y
  }
  if(drop && x$nvars == 1) {
    switch(x$vtype,
           dfcolumn = {
             return(x$df[, , drop=TRUE])
           },
           hypercolumn = {
             hc <- x$hypercolumns[[1]]
             if(x$ncases > 1) {
               hc <- as.listof(hc)
               names(hc) <- row.names(x$df)
               return(hc)
             } else {
               ha <- hc[[1]]
               return(ha)
             }
           },
           hyperatom = {
             if(x$ncases == 1) {
               # extract the hyperatom itself 
               ha <- x$hyperatoms[[1]]
               return(ha)
             } else {
               # replicate it to make a hypercolumn
               ha <- x$hyperatoms[1]
               names(ha) <- NULL
               hc <- rep(ha, x$ncases)
               hc <- as.listof(hc)
               names(hc) <- row.names(x$df)
               return(hc)
             }
           })
  }
  class(x) <- c("hyperframe", class(x))
  return(x)
}

"$.hyperframe" <- function(x,name) {
  m <- match(name, unclass(x)$vname)
  if(is.na(m))
    return(NULL)
  return(x[, name, drop=TRUE])
}

"$<-.hyperframe" <- function(x, i, value) {
  rown <- row.names(x)
  x <- as.list(x)
  dfcol <- is.atomic(value) && (is.vector(value) || is.factor(value))
  if(!dfcol && !is.null(value))
    value <- as.list(value)
  x[[i]] <- value
  y <- do.call("hyperframe", append(x, list(row.names=rown)))
  return(y)
}

"[<-.hyperframe" <- 
function (x, i, j, value)
{
  sumry <- summary(x)
  colnam <- sumry$col.names
  dimx <- sumry$dim
  die <- function(situation) {
    stop(paste("Sorry,", dQuote("[<-.hyperframe"),
               "is not yet implemented for", situation),
         call.=FALSE)
  }
  if(!missing(i))
    die("row indices")
  if(missing(j)) {
    # x[ ] <- value
    die("null indices")
  }
  if(!missing(j)) {
    # x[, j] <- value
    if(is.character(j)) {
      if(length(j) != 1)
        die("multiple columns")
      y <- get("$<-.hyperframe")(x, j, value)
      return(y)
    } else if(is.numeric(j)) {
      if(length(j) != 1 || any(j < 0))
        die("multiple columns")
      if(j <= dimx[2]) {
        jname <- colnam[j]
        y <- get("$<-.hyperframe")(x, j, value)
        return(y)
      } else if(j == dimx[2] + 1) {
        y <- cbind.hyperframe(x, value)
        return(y)
      } else 
        stop(paste("Illegal column index", j))
    } else if(is.logical(j)) {
      if(length(j) != dimx[2])
        stop(paste("Length of logical vector", paren(length(j)),
                   "does not match number of columns", paren(dimx[2])))
      if(sum(j) != 1)
        die("multiple columns")
      jname <- colnam[j]
      y <- get("$<-.hyperframe")(x, j, value)
      return(y)
    }
    else stop("Index j not understood")
  }
  return(NULL)
}

