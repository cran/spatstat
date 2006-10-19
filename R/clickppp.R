# Dominic Schuhmacher's idea

clickppp <- function(n=NULL, win=square(1), types=NULL, ..., add=FALSE, main=NULL, hook=NULL) {

  win <- as.owin(win)

  instructions <-
    if(!is.null(n)) paste("click", n, "times in window") else
    paste("add points: click left mouse button in window\n",
          "exit: click middle mouse button")
  if(is.null(main))
    main <- instructions
  
  ####  single type #########################
  if(is.null(types)) {
    plot(win, add=add, main=main)
    if(!is.null(hook))
      plot(hook, add=TRUE)
    if(!is.null(n))
      xy <- do.call("locator",
                    resolve.defaults(list(...),
                                     list(n=n, type="p")))
    else
      xy <- do.call("locator",
                    resolve.defaults(list(...),
                                     list(type="p")))
    X <- as.ppp(xy, W=win)
    Y <- X[win]
    if((ndiff <- (X$n - Y$n)) > 0)
      message(paste("deleted", ndiff,
                    ngettext(ndiff, "point", "points"),
                    "outside window"))
    return(Y)
  }
  
  ##### multitype #######################
  
  ftypes <- factor(types, levels=types)
  getem <- function(i, instr, ...) {
    main <- paste("Points of type", sQuote(i), "\n", instr)
    do.call("clickppp", resolve.defaults(list(...), list(main=main)))
  }
  # input points of type 1 
  X <- getem(ftypes[1], instructions, n=n, win=win, add=add, ..., pch=1)
  X <- X %mark% ftypes[1]
  # input points of types 2, 3, ... in turn
  for(i in 2:length(types)) {
    Xi <- getem(ftypes[i], instructions, n=n, win=win, add=add, ..., hook=X, pch=i)
    Xi <- Xi %mark% ftypes[i]
    X <- superimpose(X, Xi)
  }
  if(!add)
    plot(X, main="Final pattern")
  return(X)

}

  
