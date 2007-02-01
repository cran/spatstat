\name{nbfires}
\alias{nbfires}
\alias{nbextras}
\alias{nbw.rect}
\docType{data}
\title{Point Patterns of New Brunswick Forest Fires}
\description{
Point patterns created from yearly records, provided by the New
Brunswick Department of Natural Resources, of all fires falling
under their jurisdiction for the years 1987 to 2003 inclusive
(with the year 1988 omitted until further notice).
}
\usage{data(nbfires)}
\format{
Executing \code{data(nbfires)} creates three objects: \code{nbfires},
\code{nbextras}, and \code{nbw.rect}.

The object \code{nbfires} is a marked point pattern (an object
of class \code{"ppp"}) consisting of all of the fires in the years
1987 to 2003 inclusive, with the omission of 1988.  The marks
consist of the last two digits of the years, with the string
\code{nbfires.} prepended.  Patterns for individual years can
be extracted using the function \code{\link{split.ppp}()}.  (See
``Examples''.)

The object \code{nbextras} is a list of data frames, with
names \code{extras.87}, \ldots, \code{extras.03}.  The columns
of these data frames provide ``extra'' information about the
fires.  (See ``Details''.)

The object \code{nbw.rect} is a rectangular window which covers
central New Brunswick.  It is provided for use in illustrative and
``practice'' calculations inasmuch as the use of a rectangular
window simplifies some computations considerably.
}
\details{
The coordinates of the fire locations were provided in terms of
latitude and longitude, to the nearest minute of arc.  These were converted
to New Brunswick stereographic projection coordinates (Thomson,
Mephan and Steeves, 1977) which was the coordinate system in which
the map of New Brunswick --- which constitutes the observation
window for the pattern --- was obtained.  The conversion was done
using a \code{C} program kindly provided by Jonathan Beaudoin of
the Department of Geodesy and Geomatics, University of New Brunswick.

Finally the data and window were rescaled since the use of the
New Brunswick stereographic projection coordinate system resulted
in having to deal with coordinates which are expressed as very
large integers with a bewildering number of digits.  Amongst other
things, these huge numbers tended to create very untidy axis labels
on graphs.  The width of the bounding box of the window was made
equal to 1000 (nameless) units. In addition the lower left hand
corner of this bounding box was shifted to the origin. The height
of the bounding box was changed proportionately, resulting in a
value of approximately 959.

The window for the fire patterns comprises 6 polygonal components,
consisting of mainland New Brunswick and the 5 largest islands.
Some lakes which should form holes in the mainland component are
currently missing; this problem will be remedied in future releases.
The window was formed by ``simplifying'' the map that was originally
obtained.  The simplification consisted in reducing (using
an interactive visual technique) the number of polygon edges in
each component.  For instance the number of edges in the mainland
component was reduced from over 138,000 to 500.

For some purposes it is probably better to use a discretized
(mask type) window.  See ``Examples''.

Because of the coarseness of the coordinates of the original
data (1 minute of longitude is approximately 1 kilometer at the
latitude of New Brunswick), data entry errors, and the simplification
of the observation window, many of the original fire locations
appeared to be outside of the window.  This problem was addressed
by shifting the location of the ``outsider'' points slightly,
or deleting them, as seemed appropriate.

The columns of the data frames comprising \code{nbextras} are
\itemize{
    \item{fire.type}{A code for the type of fire with 1 =
    forest, 2 = grass, 3 = dump, and 4 = other. }
    \item{dis.date}{The discovery date of the fire, which is the
    nearest possible surrogate for the start date
    of the fire.  These columns are of class
    \code{Date} and have the format \code{yyyy-mm-dd}}
    \item{dis.time}{The discovery time (time of day) of the fire, on
    a 24 hour clock.  These
    columns are of mode character and have the format \code{HH:MM}.}
    \item{dis.julian}{The discovery date and time of the fire, expressed
    in ``Julian days'', i.e. as a decimal fraction representing the number
    of days since the beginning of the year (midnight 31 December).}
    \item{out.date}{The date on which the fire was judged to be
    ``out''.}
    \item{out.time}{The time of day on a 24 hour clock, at which
    the fire was judged to be ``out''.}
    \item{out.julian}{The date and time at which the fire was judged
    to be ``out'', expressed in Julian days.}
    \item{cause}{General cause of the fire, coded as
    1 = railroads,
    2 = unknown,
    3 = miscellaneous,
    4 = lightening,
    5 = forest industry,
    6 = incendiary,
    7 = recreation,
    8 = resident,
    9 = other industry.
    Causes 2, 4, and 6 are designated as ``final'' by the New
    Brunswick Department of Natural Resources, meaning
    (it seems) ``that's all there is to it''.  Other causes
    are apparently intended to be refined by being combined
    with ``source of ignition''.}
    \item{ign.src}{Source of ignition, coded as
    1 = cigarette/match/pipe/ashes,
    2 = burning without a permit,
    3 = burning with a permit,
    4 = prescribed burn,
    5 = wood spark,
    6 = machine spark,
    7 = campfire,
    8 = chainsaw,
    9 = machinery,
    10 = vehicle accident,
    11 = railroad accident,
    12 = wheelbox on railcars,
    13 = hot flakes off railcar wheels,
    14 = dump (i.e. fire escaping from a dump),
    15 = ashes (briquettes, burning garbage, etc.)}
    \item{fnl.size}{The final size of the fire (area burned)
    in hectares, to the nearest 10th hectare.}
}
Note that due to data entry errors some of the ``out dates'' and
``out times'' in the original data sets were actually \emph{earlier}
than the corresponding ``discovery dates'' and ``discover times''.
In such cases all corresponding entries of the \code{extras} data
frame (i.e. \code{dis.date}, \code{dis.time}, \code{dis.julian},
\code{out.date}, \code{out.time}, and \code{out.julian}) were set
equal to \code{NA}.  Also, some of the dates and times were missing
(equal to \code{NA}) in the original data sets.
}
\source{
  The data were kindly provided by the New Brunswick Department
  of Natural Resources.  Special thanks are due to Jefferey Betts for
  a great deal of assistance.
}
\references{
Turner, Rolf (2007)
Point patterns of forest fire locations.
To appear in \emph{Environmental and Ecological Statistics}.

Thomson, D. B., Mephan, M. P., and Steeves, R. R. (1977)
The stereographic double projection.
Technical Report 46, University of New Brunswick,
Fredericton, N. B., Canada
URL: \code{gge.unb.ca/Pubs/Pubs.html}.
}
\examples{
\dontrun{
data(nbfires)
X <- split(nbfires) # Create a list of yearly point patterns.
attach(X)           # Make the individual point patterns accessible.
attach(nbextras)    # Make the individual `extras' accessible.
ftyp <- extras.00$fire.type
Y.00 <- nbfires.00[ftyp==1 | ftyp==2] # Pick out forest and grass fires.
stt  <- extras.00$dis.julian[ftyp==1 | ftyp==2] 
fin  <- extras.00$out.julian[ftyp==1 | ftyp==2] 
Y.00 <- setmarks(Y.00,fin-stt) # Mark the pattern with fire duration.
plot(Y.00)
#
nbw.mask <- as.mask(nbfires$window, dimyx=500)
plot(nbw.mask, col=c("green", "white"))
plot(nbfires$window, border="red", add=TRUE)
#
plot(unmark(Y.00)[nbw.rect], add=TRUE)
plot(nbw.rect,add=TRUE,border="blue")
#
K.00 <- Kest(Y.00)
plot(K.00)
}
}
\keyword{datasets}
\keyword{spatial}
