#
#     options.R
#
#     Spatstat Options
#
#    $Revision: 1.28 $   $Date: 2010/06/05 10:50:27 $
#
#

".Spatstat.Options" <- list()

reset.spatstat.options <- function() {
  .Spatstat.Options <<- lapply(.Spat.Stat.Opt.Table,
                               function(z) { z$default })
  invisible(.Spatstat.Options)  
}

".Spat.Stat.Opt.Table" <-
  list(
       npixel=list(
         default=100,
         check=function(x){
           is.numeric(x) && (length(x) %in% c(1,2)) && is.finite(x) &&
           all(x == ceiling(x)) && all(x > 1) 
         },
         valid="an integer, or a pair of integers, greater than 1"
        ),
       maxedgewt=list(
         default=100,
         check=function(x){
           is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1
         },
         valid="a finite numeric value, not less than 1"
       ),
       par.binary=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.persp=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.points=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       par.contour=list(
         default=list(),
         check=is.list,
         valid="a list"
         ),
       image.colfun=list(
         default=function(n){topo.colors(n)},
         check=function(x) {
           is.function(x) && length(formals(x)) > 0 && all(is.character(x(42)))
         },
         valid="a function that returns character values"
         ),
       ndummy.min=list(
         default=25,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1
         },
         valid="a single integer, greater than 1"
       ),
       dupC = list(
         default=FALSE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       progress = list(
         default="tty",
         check=function(x){ x %in% c("tty", "txtbar") },
         valid=paste("one of the strings", dQuote("tty"),
           "or", dQuote("txtbar"))
         ),
       checkpolygons = list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       checksegments = list(
         default=TRUE,
         check=function(x) { is.logical(x) && length(x) == 1},
         valid="a single logical value"
         ),
       ngrid.disc=list(
         default=128,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1
         },
         valid="a single integer, greater than 1"
       ),
       gpclib=list(
         default=FALSE,
         check=function(x) {
           if(!(is.logical(x) && length(x) == 1))
             return(FALSE)
           if(x && !require(gpclib)) {
             warning("Cannot set gpclib=TRUE: package gpclib is not installed")
             return(FALSE)
           }
           return(TRUE)
         },
         valid="a single logical value"
         ),
       maxmatrix=list(
         default=2^24, # 16,777,216
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
       ),
       huge.npoints=list(
         default=1e6,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && (x == ceiling(x)) && x > 1024
         },
         valid="a single integer, greater than 1024"
       ),
       expand=list(
         default=2,
         check=function(x) {
           is.numeric(x) && length(x) == 1 && x > 1
         },
         valid="a single numeric value, greater than 1"
       )
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
# validate new values
    for(i in seq(assignto)) {
      nama <- assignto[i]
      valo <- called[[i]]
      entry <- .Spat.Stat.Opt.Table[[nama]]
      ok <- do.call(entry$check, list(valo))
      if(!ok)
        stop(paste("Parameter", dQuote(nama), "should be",
                   entry$valid))
    }
# reassign
    changed <- .Spatstat.Options[assignto]
    .Spatstat.Options[assignto] <<- called
# return 
    invisible(changed)
}

