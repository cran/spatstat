#
#	interact.S
#
#
#	$Revision: 1.7 $	$Date: 2004/06/08 10:48:08 $
#
#	Class 'interact' representing the interpoint interaction
#               of a point process model
#              (e.g. Strauss process with a given threshold r)
#
#       Class 'isf' representing a generic interaction structure
#              (e.g. pairwise interactions)
#
#	These do NOT specify the "trend" part of the model,
#	only the "interaction" component.
#
#               The analogy is:
#
#                       glm()             ppm()
#
#                       model formula     trend formula
#
#                       family            interaction
#
#               That is, the 'systematic' trend part of a point process
#               model is specified by a 'trend' formula argument to ppm(),
#               and the interpoint interaction is specified as an 'interact'
#               object.
#
#       You only need to know about these classes if you want to
#       implement a new point process model.
#
#       THE DISTINCTION:
#       An object of class 'isf' describes an interaction structure
#       e.g. pairwise interaction, triple interaction,
#       pairwise-with-saturation, Dirichlet interaction.
#       Think of it as determining the "order" of interaction
#       but not the specific interaction potential function.
#
#       An object of class 'interact' completely defines the interpoint
#       interactions in a specific point process model, except for the
#       regular parameters of the interaction, which are to be estimated
#       by ppm() or otherwise. An 'interact' object specifies the values
#       of all the 'nuisance' or 'irregular' parameters. An example
#       is the Strauss process with a given, fixed threshold r
#       but with the parameters beta and gamma undetermined.
#
#       DETAILS:
#
#       An object of class 'isf' contains the following:
#
#	     $name               Name of the interaction structure         
#                                        e.g. "pairwise"
#
#	     $print		 How to 'print()' this object
#				 [A function; invoked by the 'print' method
#                                 'print.isf()']
#
#            $eval               A function which evaluates the canonical
#                                sufficient statistic for an interaction
#                                of this general class (e.g. any pairwise
#                                interaction.)
#
#       If lambda(u,X) denotes the conditional intensity at a point u
#       for the point pattern X, then we assume
#                  log lambda(u, X) = theta . S(u,X)
#       where theta is the vector of regular parameters,
#       and we call S(u,X) the sufficient statistic.
#
#       A typical calling sequence for the $eval function is
#
#            (f$eval)(X, U, E, potentials, potargs, correction)
#
#       where X is the data point pattern, U is the list of points u
#       at which the sufficient statistic S(u,X) is to be evaluated,
#       E is a logical matrix equivalent to (X[i] == U[j]),
#       $potentials defines the specific potential function(s) and
#       $potargs contains any nuisance/irregular parameters of these
#       potentials [the $potargs are passed to the $potentials without
#       needing to be understood by $eval.]
#       $correction is the name of the edge correction method.
#
#
#       An object of class 'interact' contains the following:
#
#
#            $name               Name of the specific potential
#                                        e.g. "Strauss"
#
#            $family              Object of class "isf" describing
#                                the interaction structure
#
#            $pot	         The interaction potential function(s)
#                                -- usually a function or list of functions.
#                                (passed as an argument to $family$eval)
#
#            $par                list of any nuisance/irregular parameters
#                                (passed as an argument to $family$eval)
#
#            $parnames           vector of long names/descriptions
#                                of the parameters in 'par'
#
#            $init()             initialisation action
#                                or NULL indicating none required
#
#            $update()           A function to modify $par
#                                [Invoked by 'update.interact()']
#                                or NULL indicating a default action
#
#	     $print		 How to 'print()' this object
#				 [Invoked by 'print' method 'print.interact()']
#                                or NULL indicating a default action
#
# --------------------------------------------------------------------------

print.isf <- function(x, ...) {
  verifyclass(x, "isf")
  if(!is.null(x$print))
    (x$print)(x)
  invisible(NULL)
}

print.interact <- function(x, ...) {
  verifyclass(x, "interact")
  if(!is.null(x$print))
     (x$print)(x)
  else {
    # default
    print.isf(x$family)
    cat(paste("Interaction:", x$name, "\n"))
    # just print the parameter names and their values
    cat(paste(x$parnames, ":\t", x$par, "\n", sep=""))
  }
  invisible(NULL)
}

update.interact <- function(object, ...) {
  verifyclass(object, "interact")
  if(!is.null(object$update))
    (object$update)(object, ...)
  else {
    # Default
    # just match the arguments in "..."
    # with those in object$par and update them
    want <- list(...)
    m <- match(names(want),names(object$par))
    nbg <- is.na(m)
    if(any(nbg)) {
        which <- paste((names(want))[nbg])
        warning(paste("Arguments not matched: ", which))
    }
    m <- m[!nbg]
    object$par[m] <- u
    # call object's own initialisation routine
    if(!is.null(object$init))
      (object$init)(object)
    object
  }    
}

