\name{spatstat-internal}
\alias{affinexy}
\alias{apply23sum}
\alias{area.xypolygon}       
\alias{as.breakpts}
\alias{as.polygonal}
\alias{bdrylength.xypolygon} 
\alias{breakpts}
\alias{breakpts.from.r}
\alias{checkfields}          
\alias{check.hist.lengths}
\alias{countingweights}
\alias{crossdist}            
\alias{default.ngrid}
\alias{discretise}
\alias{distpl}               
\alias{distppl}
\alias{distppll}
\alias{equals.quad}          
\alias{erode.mask}
\alias{even.breaks.owin}
\alias{exactdt}              
\alias{exactPdt}
\alias{getfields}
\alias{gridindex}            
\alias{handle.r.b.args}
\alias{identical.formulae}
\alias{inside.xypolygon}     
\alias{is.data}
\alias{is.marked}
\alias{is.marked.default}    
\alias{is.marked.ppp}
\alias{is.poisson.ppm}
\alias{is.stationary.ppm}    
\alias{km.rs}
\alias{Kount}
\alias{Kwtsum}               
\alias{make.even.breaks}
\alias{markspace.integral}
\alias{marks.quad}           
\alias{matcolall}
\alias{matcolany}
\alias{matcolsum}            
\alias{matrowall}
\alias{matrowany}
\alias{matrowsum}            
\alias{MultiPair.checkmatrix}
\alias{multipair.family}
\alias{nndist}
\alias{no.trend.ppm}
\alias{n.quad}               
\alias{overlap.owin}
\alias{overlap.trapezium}
\alias{overlap.xypolygon}    
\alias{pairdist}
\alias{polynom}
\alias{print.interact}       
\alias{print.isf}
\alias{quad}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{replicate}
\alias{reverse.xypolygon}    
\alias{rotxy}
\alias{runifdisc}
\alias{runifpoispp}          
\alias{runifrect}
\alias{scanpp}
\alias{shiftxy}              
\alias{sympoly}
\alias{termsinformula}
\alias{tilecentroids}        
\alias{tilemiddles}
\alias{try.FGJKest}
\alias{try.Kcross}           
\alias{unit.square}
\alias{update.interact}
\alias{validate.mask}        
\alias{verifyclass}
\alias{verify.xypolygon}
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
as.polygonal(W)
bdrylength.xypolygon(polly)
breakpts(val, maxi, even = F, npos = NULL, step = NULL)
breakpts.from.r(r)
checkfields(X,L)          
check.hist.lengths(hist,breaks)
countingweights(id, areas, check = T)
crossdist(x1, y1, x2, y2)            
default.ngrid(X)
discretise(x, xrange, nx)
distpl(p, l)               
distppl(p, l)
distppll(p, l)
equals.quad(Q)          
erode.mask(w,r)
even.breaks.owin(w)
exactdt(X, \dots)              
exactPdt(im)
getfields(X, L, fatal = T)
gridindex(x, y, xrange, yrange, nx, ny)            
handle.r.b.args(r = NULL, breaks = NULL, window, eps = NULL)
identical.formulae(a,b)
inside.xypolygon(pts, polly, test01 = T)     
is.data(Q)
is.marked(X, \dots)
is.marked.default(\dots)  
is.marked.ppp(X, na.action="warn")
is.poisson.ppm(x)
is.stationary.ppm(x)    
km.rs(o, cc, d, breaks)
Kount(d, b, breaks, slow = F)
Kwtsum(d, b, wi, wj, breaks)               
make.even.breaks(bmax, npos, bstep)
markspace.integral(X)
marks.quad(Q)           
matcolall(x)
matcolany(x)
matcolsum(x)            
matrowall(x)
matrowany(x)
matrowsum(x)            
MultiPair.checkmatrix(mat, n, name)
multipair.family
nndist(x, y)
no.trend.ppm(x)
n.quad(Q)               
overlap.owin(A,B)
overlap.trapezium(xa, ya, xb, yb, verb = F)
overlap.xypolygon(P, Q)   
pairdist(x, y)
polynom(x, \dots)
print.interact(x, \dots)       
print.isf(x, \dots)
quad(data, dummy, w)
quadscheme.replicated(data, dummy = default.dummy(data), method = "grid", \dots)
quadscheme.spatial(data, dummy = default.dummy(data), method = "grid", \dots)
replicate(pp, markset, fac = T)
reverse.xypolygon(p)    
rotxy(X, angle = pi/2)
runifdisc(n, r = 1, x = 0, y = 0)
runifpoispp(lambda, win = owin(c(0, 1), c(0, 1)))
runifrect(n, win = owin(c(0, 1), c(0, 1)))
scanpp(filename, window, header = T, dir = "", multitype = F)
shiftxy(X, vec = c(0, 0))              
sympoly(x, y, n)
termsinformula(x)
tilecentroids(W, nx, ny)        
tilemiddles(W, nx, ny)
try.FGJKest(niter = 20, lambda = 25, r = seq(0, sqrt(2), 0.02), eps = 0.01, slow = F)
try.Kcross(niter = 20, lambda1 = 25, lambda2 = 25, r = seq(0, 1, 0.02), R = 0.2)           
unit.square()
update.interact(object, \dots)
validate.mask(w, fatal=T)        
verifyclass(X, C, N = deparse(substitute(X)), fatal = T)
verify.xypolygon(p, fatal=T)
w.quad(Q)               
x.quad(Q)
y.quad(Q)               
}

\details{
  These are not to be called by the user.
}
\keyword{internal}
