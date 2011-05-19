#
# split.ppp.R
#
# $Revision: 1.15 $ $Date: 2011/05/18 09:14:29 $
#
# split.ppp and "split<-.ppp"
#
#########################################

split.ppp <- function(x, f = marks(x), drop=FALSE, un=NULL, ...) {
  verifyclass(x, "ppp")
  mf <- markformat(x)
  
  if(is.null(un))
    un <- missing(f) && (mf != "dataframe")

  if(missing(f)) {
    # f defaults to marks of x
    switch(mf,
           none={
             stop("f is missing and there are no marks")
           },
           vector={
             if(!is.multitype(x)) 
               stop("f is missing and the pattern is not multitype")
             f <- fsplit <- marks(x)
           },
           dataframe={
             f <- fsplit <- firstfactor(marks(x))
             if(is.null(f))
               stop("Data frame of marks contains no factors")
           })
  } else{
    # f was given
    fsplit <- f
    if(inherits(f, "im")) {
      # f is an image: determine the grouping
      fsplit <- tess(image=f)
      f <- marks(cut(x, fsplit))
    } else if(inherits(f, "tess")) {
      # f is a tessellation: determine the grouping
      f <- marks(cut(x, fsplit))
    } else if(is.character(f) && length(f) == 1) {
      # f is the name of a column of marks
      marx <- marks(x)
      if(is.data.frame(marx) && (f %in% names(marx))) 
        fsplit <- f <- marx[[f]]
      else
        stop(paste("The name", sQuote(f), "does not match any column of marks"))
    }
    # validate
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  }

  lev <- levels(f)
  if(drop) {
    lev <- lev[table(f) > 0]
    fsplit <- fsplit[f %in% lev]
  }

  # split the data
  out <- list()
  for(l in lev) 
    out[[paste(l)]] <- x[f == l]
  
  if(un)
     out <- lapply(out, unmark)
  if(inherits(fsplit, "tess")) {
    til <- tiles(fsplit)
    for(i in seq_along(out))
      out[[i]]$window <- til[[i]]
  }
  class(out) <- c("splitppp", class(out))
  attr(out, "fsplit") <- fsplit
  return(out)
}

"split<-.ppp" <- function(x, f=marks(x), drop=FALSE, un=missing(f), 
                          ..., value) {
  verifyclass(x, "ppp")
  stopifnot(is.list(value))
  if(!all(unlist(lapply(value, is.ppp))))
    stop(paste("Each entry of", sQuote("value"),
               "must be a point pattern"))

  ismark <- unlist(lapply(value, is.marked))
  if(any(ismark) && !all(ismark))
    stop(paste("Some entries of",
               sQuote("value"),
               "are marked, and others are unmarked"))
  vmarked <- all(ismark)

  # evaluate `un' before assigning value of 'f'
  un <- un
  
  if(!missing(f)) {
    fsplit <- f
    if(inherits(f, "tess"))
      f <- marks(cut(x, f))
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } else {
    if(is.multitype(x))
      f <- fsplit <- marks(x)
    else
      stop("f is missing and there is no sensible default")
  }
  lev <- levels(f)
  if(!drop) 
    levtype <- "levels of f"
  else {
    levtype <- "levels which f actually takes"
    lev <- lev[table(f) > 0]
    fsplit <- fsplit[f %in% lev]
  }
  if(length(value) != length(lev))
      stop(paste("length of", sQuote("value"),
                 "should equal the number of",
                 levtype))

  # ensure value[[i]] is associated with lev[i]
  if(!is.null(names(value))) {
    if(!all(names(value) %in% paste(levels(f))))
      stop(paste("names of", sQuote("value"), "should be levels of f"))
    value <- value[lev]
  }
  names(value) <- NULL
  
  # restore the marks, if they were discarded
  if(un && is.marked(x)) {
    if(vmarked)
      warning(paste(sQuote("value"), "contains marked point patterns:",
                    "this is inconsistent with un=TRUE; marks ignored."))
    for(i in seq_along(value)) 
      value[[i]] <- value[[i]] %mark% factor(lev[i], levels=levels(f))
  }

  #out <- superimpose(value, W=x$window)
  out <- do.call(superimpose,c(value,list(W=x$window)))
  return(out)
}


print.splitppp <- function(x, ...) {
  f <- attr(x, "fsplit")
  cat(paste("Point pattern split by",
            if(inherits(f, "tess")) "tessellation" else "factor",
            "\n"))
  nam <- names(x)
  for(i in seq_along(x)) {
    cat(paste("\n", nam[i], ":\n", sep=""))
    print(x[[i]])
  }
  return(invisible(NULL))
}

summary.splitppp <- function(object, ...) {
  x <- lapply(object, summary, ...)
  class(x) <- "summary.splitppp"
  x
}

print.summary.splitppp <- function(x, ...) {
  class(x) <- "listof"
  print(x)
  invisible(NULL)
}

"[.splitppp" <- function(x, ...) {
  f <- attr(x, "fsplit")
  # invoke list method on x
  class(x) <- "list"
  y <- x[...]
  # then make it a 'splitppp' object too
  class(y) <- c("splitppp", class(y))
  attr(y, "fsplit") <- f[...]
  y
}

"[<-.splitppp" <- function(x, ..., value) {
  f <- attr(x, "fsplit")
  # invoke list method
  class(x) <- "list"
  x[...] <- value
  # then make it a 'splitppp' object too
  class(x) <- c("splitppp", class(x))
  attr(x, "fsplit") <- f
  x
}
  
