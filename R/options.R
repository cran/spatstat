#
#     options.R
#
#     Spatstat Options
#
#    $Revision: 1.8 $   $Date: 2008/04/02 13:32:43 $
#
#
".Spatstat.Options" <-
  list(npixel = 100,
       maxedgewt=100.0,
       par.binary=list(),
       par.persp=list(),
       par.points=list(),
       par.contour=list(),
       ndummy.min=25,
       dupC = FALSE
       )

"spatstat.options" <-
function (...) 
{
    called <- list(...)    

    if(length(called) == 0)
    	return(.Spatstat.Options)

    if(is.null(names(called)) && length(called)==1) {
      # spatstat.options(x) 
      x <- called[[1]]
      if(is.null(x))
        return(.Spatstat.Options)  # spatstat.options(NULL)
      if(is.list(x))
        called <- x 
    }
    
    if(is.null(names(called))) {
        # spatstat.options("par1", "par2", ...)
	ischar <- unlist(lapply(called, is.character))
	if(all(ischar)) {
		choices <- unlist(called)
		ok <- choices %in% names(.Spatstat.Options)
		if(!all(ok))
                  stop(paste("Unrecognised option(s):", called[!ok]))
                if(length(called) == 1)
                  return(.Spatstat.Options[[choices]])
                else
                  return(.Spatstat.Options[choices])
	} else {
	   wrong <- called[!ischar]
	   offending <- unlist(lapply(wrong,
	   		function(x) { y <- x;
	     		deparse(substitute(y)) }))
	   offending <- paste(offending, collapse=",")
           stop(paste("Unrecognised mode of argument(s) [",
		offending,
	   "]: should be character string or name=value pair"))
    	}
    }
# spatstat.options(name=value, name2=value2,...)
    assignto <- names(called)
    if (is.null(assignto) || any(assignto == "")) 
        stop("options must all be identified by name=value")
    ok <- assignto %in% names(.Spatstat.Options)
    if(!all(ok))
	stop(paste("Unrecognised option(s):", assignto[!ok]))
# reassign
    changed <- .Spatstat.Options[assignto]
    .Spatstat.Options[assignto] <<- called
# return 
    invisible(changed)
}

