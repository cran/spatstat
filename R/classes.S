#
#
#	classes.S
#
#	$Revision: 1.3 $	$Date: 2000/07/11 10:53:19 $
#
#	Generic utilities for classes
#
#	These are supposed to work under BOTH Splus3.4 and Splus5.1
#
#
#--------------------------------------------------------------------------

verifyclass <- function(X, C, N=deparse(substitute(X)), fatal=T) {
  if(is.null(class(X)) || class(X) != C) {
    if(fatal) {
        gripe <- paste("argument \'", N,
                       "\' is not of class \'", C, "\'", sep="")
	stop(gripe)
    } else 
	return(F)
  }
  return(T)
}

#--------------------------------------------------------------------------

checkfields <- function(X, L) {
	  # X is a list, L is a vector of strings
	  # Checks for presence of field named L[i] for all i
	return(all(!is.na(match(L,names(X)))))
}

getfields <- function(X, L, fatal=T) {
	  # X is a list, L is a vector of strings
	  # Extracts all fields with names L[i] from list X
	  # Checks for presence of all desired fields
	  # Returns the sublist of X with fields named L[i]
	absent <- is.na(match(L, names(X)))
	if(any(absent)) {
		gripe <- paste("Needed the following components:",
				paste(L, collapse=", "),
				"\nThese ones were missing: ",
				paste(L[absent], collapse=", "))
		if(fatal)
			stop(gripe)
		else 
			warning(gripe)
	} 
	return(X[L[!absent]])
}



