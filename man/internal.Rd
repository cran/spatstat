\name{spatstat-internal}
\alias{affinexy}
\alias{apply23sum}
\alias{area.xypolygon}       
\alias{as.breakpts}
\alias{as.fv}
\alias{as.polygonal}
\alias{bind.fv}
\alias{bdrylength.xypolygon} 
\alias{breakpts}
\alias{breakpts.from.r}
\alias{cartesian}
\alias{checkfields}          
\alias{check.hist.lengths}
\alias{conspire}
\alias{contour.im}
\alias{countingweights}
\alias{crossdist}            
\alias{default.ngrid}
\alias{default.ntile}
\alias{discretise}
\alias{distpl}               
\alias{distppl}
\alias{distppll}
\alias{edge.Ripley}
\alias{edge.Trans}
\alias{equals.quad}          
\alias{erode.mask}
\alias{even.breaks.owin}
\alias{exactdt}              
\alias{exactPdt}
\alias{fasp}
\alias{fv}
\alias{getfields}
\alias{gridindex}            
\alias{handle.r.b.args}
\alias{identical.formulae}
\alias{image.im}            
\alias{inside.xypolygon}     
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}     
\alias{is.data}
\alias{is.fv}
\alias{is.marked.default}    
\alias{is.poisson.ppm}
\alias{is.stationary.ppm}    
\alias{Kount}
\alias{Kwtsum}               
\alias{lookup.im}
\alias{make.even.breaks}
\alias{markspace.integral}
\alias{marks.quad}           
\alias{matcolall}
\alias{matcolany}
\alias{matcolsum}            
\alias{matrowall}
\alias{matrowany}
\alias{matrowsum}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{mkcor}
\alias{mpl}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{MultiPair.checkmatrix}
\alias{multipair.family}
\alias{nearest.pixel}
\alias{no.trend.ppm}
\alias{n.quad}               
\alias{overlap.owin}
\alias{overlap.trapezium}
\alias{overlap.xypolygon}    
\alias{param.quad}
\alias{persp.im}
\alias{pixellate}
\alias{plot.im}
\alias{polynom}
\alias{print.fasp}       
\alias{print.fv}       
\alias{print.interact}       
\alias{print.isf}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{quad}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{rasterx.im}
\alias{rastery.im}
\alias{reverse.xypolygon}    
\alias{rotxy}
\alias{runifdisc}
\alias{runifpoispp}          
\alias{runifrect}
\alias{rmhmodel.ppm}
\alias{rmpoint.I.allim}
\alias{rpoint.multi}
\alias{scanpp}
\alias{second.moment.calc}
\alias{setcov}
\alias{shiftxy}              
\alias{spatstat.rawdata.location}
\alias{sympoly}
\alias{termsinformula}
\alias{tilecentroids}        
\alias{tilemiddles}
\alias{try.FGJKest}
\alias{try.Kcross}           
\alias{update.interact}
\alias{validate.mask}        
\alias{verifyclass}
\alias{verify.xypolygon}
\alias{whist}               
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}               
\title{Internal spatstat functions}
\description{
  Internal spatstat functions.
}
\usage{
affinexy(X, mat = diag(c(1, 1)), vec = c(0, 0))
apply23sum(x)
area.xypolygon(polly)       
as.breakpts(\dots)
as.fv(x)
as.polygonal(W)
bind.fv(x, y, labl, desc, preferred)
bdrylength.xypolygon(polly)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
cartesian(pp, markset, fac = TRUE)
checkfields(X,L)          
check.hist.lengths(hist,breaks)
conspire(\dots)
contour.im(x, \dots, xlab="x", ylab="y")
countingweights(id, areas, check = TRUE)
crossdist(x1, y1, x2, y2)            
default.ngrid(X)
default.ntile(X)
discretise(x, xrange, nx)
distpl(p, l)               
distppl(p, l)
distppll(p, l)
edge.Ripley(X, r, W)
edge.Trans(X, Y, W, exact, trim)
equals.quad(Q)          
erode.mask(w,r)
even.breaks.owin(w)
exactdt(X, \dots)              
exactPdt(im)
fv(x, argu, ylab, valu, fmla, alim, labl, desc)
fasp((fns, titles, formulae, which, dataname, title)
getfields(X, L, fatal = TRUE)
gridindex(x, y, xrange, yrange, nx, ny)            
handle.r.b.args(r = NULL, breaks = NULL, window, eps = NULL)
identical.formulae(x,y)
image.im(x, \dots, xlab="x", ylab="y")
inside.xypolygon(pts, polly, test01 = TRUE)     
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
is.data(Q)
is.fv(x)
is.marked(X, \dots)
is.marked.default(\dots)  
is.marked.ppp(X, na.action="warn", \dots)
is.poisson.ppm(x)
is.stationary.ppm(x)    
km.rs(o, cc, d, breaks)
Kount(d, b, breaks, slow = FALSE)
Kwtsum(d, b, weight, breaks, slow)
lookup.im(im, x, y, naok)
make.even.breaks(bmax, npos, bstep)
markspace.integral(X)
marks.quad(Q)           
matcolall(x)
matcolany(x)
matcolsum(x)            
matrowall(x)
matrowany(x)
matrowsum(x)
meanX.owin(w)            
meanY.owin(w)
mkcor(d, ff, wt, Ef, rvals, method="smrep", ..., nwtsteps=500)
mpl(Q, trend, interaction, data, correction, rbord, use.gam)
mpl.engine(Q, trend, interaction, covariates, correction,
	 rbord, use.gam)
mpl.get.covariates(covariates, locations, type)
mpl.prepare(Q, X, P, trend, interaction, covariates,
            want.trend, want.inter, correction, rbord)
MultiPair.checkmatrix(mat, n, name)
multipair.family
nearest.pixel(x, y, im)
no.trend.ppm(x)
n.quad(Q)               
overlap.owin(A,B)
overlap.trapezium(xa, ya, xb, yb, verb = FALSE)
overlap.xypolygon(P, Q)   
param.quad(Q)
persp.im(x, \dots, xlab="x", ylab="y")
pixellate(x, \dots, weights)
plot.im(x, \dots, xlab="x", ylab="y")
polynom(x, \dots)
print.fasp(x, \dots)
print.fv(x, \dots)
print.interact(x, \dots)       
print.isf(x, \dots)
print.summary.owin(x, \dots)
print.summary.ppp(x, \dots, dp=3)
quad(data, dummy, w, param)
quadscheme.replicated(data, dummy, method = "grid", \dots)
quadscheme.spatial(data, dummy, method = "grid", \dots)
rasterx.im(x)
rastery.im(x)
reverse.xypolygon(p)    
rotxy(X, angle = pi/2)
rmhmodel.ppm(model, win, \dots, verbose, project)
rmpoint.I.allim(n, f, types)
rpoint.multi(n, f, fmax=NULL, marks = NULL, win = unit.square(), giveup = 1000, verbose = FALSE)
runifdisc(n, r = 1, x = 0, y = 0)
runifpoispp(lambda, win = owin(c(0, 1), c(0, 1)))
runifrect(n, win = owin(c(0, 1), c(0, 1)))
scanpp(filename, window, header = TRUE, dir = "", multitype = FALSE)
second.moment.calc(x, sigma, edge=TRUE, what="Kmeasure", debug=FALSE, \dots)
setcov(W)
shiftxy(X, vec = c(0, 0))
spatstat.rawdata.location()             
sympoly(x, y, n)
termsinformula(x)
tilecentroids(W, nx, ny)        
tilemiddles(W, nx, ny)
try.FGJKest(niter = 20, lambda = 25, r = seq(0, sqrt(2), 0.02), eps = 0.01, slow = FALSE)
try.Kcross(niter = 20, lambda1 = 25, lambda2 = 25, r = seq(0, 1, 0.02), R = 0.2)           
update.interact(object, \dots)
validate.mask(w, fatal=TRUE)        
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
verify.xypolygon(p, fatal=TRUE)
whist(x,breaks,weights)
w.quad(Q)               
x.quad(Q)
y.quad(Q)               
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
