\name{spatstat-internal}
\alias{affinexy}
\alias{apply23sum}
\alias{area.xypolygon}       
\alias{AsymmHausdorff.psp}
\alias{as.breakpts}
\alias{as.fv}
\alias{as.polygonal}
\alias{beachcolours}
\alias{bind.fv}
\alias{bdrylength.xypolygon} 
\alias{breakpts}
\alias{breakpts.from.r}
\alias{cartesian}
\alias{cellmiddles}
\alias{checkfields}          
\alias{check.hist.lengths}
\alias{check.named.vector}
\alias{clip.psp}
\alias{cliprect.psp}
\alias{clippoly.psp}
\alias{closepairs}
\alias{crosspairs}
\alias{cobble.xy}
\alias{commasep}
\alias{contour.im}
\alias{countingweights}
\alias{damaged.ppm}
\alias{default.expand}
\alias{default.clipwindow}
\alias{default.n.tiling}
\alias{default.ntile}
\alias{densityhack}
\alias{diagnose.ppm.engine}
\alias{discretise}
\alias{distpl}               
\alias{distppl}
\alias{distppll}
\alias{divisors}
\alias{do.call.matched}
\alias{edge.Ripley}
\alias{edge.Trans}
\alias{ensure2vector}
\alias{equals.quad}          
\alias{eratosthenes}
\alias{erode.mask}
\alias{even.breaks.owin}
\alias{exactdt}              
\alias{exactPdt}
\alias{fasp}
\alias{FormatFaspFormulae}
\alias{fv}
\alias{greatest.common.divisor}
\alias{getfields}
\alias{getglmfit}
\alias{gridindex}            
\alias{handle.r.b.args}
\alias{handle.rshift.args}
\alias{ho.engine}
\alias{identical.formulae}
\alias{inside.xypolygon}     
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}     
\alias{is.cadlag}
\alias{is.data}
\alias{is.fv}
\alias{is.marked.default}    
\alias{is.prime}
\alias{is.poisson.ppm}
\alias{is.stationary.ppm}    
\alias{Kborder.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{Kmulti.inhom}               
\alias{killinteraction}
\alias{least.common.multiple}
\alias{levels<-.im}
\alias{lookup.im}
\alias{make.even.breaks}
\alias{markspace.integral}
\alias{marks.quad}           
\alias{matcolall}
\alias{matcolany}
\alias{matcolsum}
\alias{matrixsample}
\alias{matrowall}
\alias{matrowany}
\alias{matrowsum}
\alias{maxflow}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{mkcor}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{MultiPair.checkmatrix}
\alias{nearest.pixel}
\alias{no.trend.ppm}
\alias{n.quad}               
\alias{overlap.owin}
\alias{overlap.trapezium}
\alias{overlap.xypolygon}    
\alias{owinpolycheck}
\alias{param.quad}
\alias{pixellate}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.diagppm}
\alias{plot.minconfit}
\alias{plot.pppmatching}
\alias{plot.profilepl}
\alias{plot.qqppm}
\alias{plot.quadratcount}
\alias{plot.quadrattest}
\alias{polynom}
\alias{pppdist.prohorov}
\alias{pppmatching}
\alias{primefactors}       
\alias{print.diagppm}
\alias{print.envelope}
\alias{print.fasp}       
\alias{print.fv}       
\alias{print.interact}       
\alias{print.isf}
\alias{print.minconfit}
\alias{print.plotppm}
\alias{print.pppmatching}
\alias{print.profilepl}
\alias{print.qqppm}
\alias{print.rmhcontrol}
\alias{print.rmhmodel}
\alias{print.rmhstart}
\alias{print.rmhseed}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{print.summary.psp}
\alias{quad}
\alias{quadrat.breaks}
\alias{quadrat.count.engine}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{rasterx.im}
\alias{rastery.im}
\alias{rebound.ppp}
\alias{rebound.psp}
\alias{rebound.owin}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{resolve.defaults}
\alias{reverse.xypolygon}    
\alias{rmax.rule}
\alias{rotxy}
\alias{runifdisc}
\alias{runifpoispp}          
\alias{runifrect}
\alias{rmhResolveTypes}
\alias{rmhcontrol.rmhcontrol}
\alias{rmhcontrol.list}
\alias{rmhEngine}
\alias{rmhmodel.rmhmodel}
\alias{rmhmodel.list}
\alias{rmhseed}
\alias{rmhstart.rmhstart}
\alias{rmhstart.list}
\alias{rmhmodel.ppm}
\alias{rmpoint.I.allim}
\alias{rpoint.multi}
\alias{second.moment.calc}
\alias{shiftxy}              
\alias{spatstat.rawdata.location}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{stratrand}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summary.envelope}
\alias{summary.pppmatching}
\alias{sympoly}
\alias{termsinformula}
\alias{tilecentroids}        
\alias{trim.rectangle}
\alias{update.interact}
\alias{validate.mask}        
\alias{validate.quad}        
\alias{variablesinformula}        
\alias{verifyclass}
\alias{verify.xypolygon}
\alias{whist}               
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}
\alias{xypolygon2psp}
\alias{xypolyselfint}
\title{Internal spatstat functions}
\description{
  Internal spatstat functions.
}
\usage{
affinexy(X, mat = diag(c(1, 1)), vec = c(0, 0))
apply23sum(x)
area.xypolygon(polly)
AsymmHausdorff.psp(X, Y, method="Fortran")
as.breakpts(\dots)
as.fv(x)
as.polygonal(W)
beachcolours(heightrange, sealevel, monochrome, ncolours)
bind.fv(x, y, labl, desc, preferred)
bdrylength.xypolygon(polly)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
cartesian(pp, markset, fac = TRUE)
cellmiddles(W, nx, ny, npix, gi)
checkfields(X,L)          
check.hist.lengths(hist,breaks)
check.named.vector(x, nam)
clip.psp(x, window, check=TRUE)
cliprect.psp(x, window)
clippoly.psp(s, window)
closepairs(X,rmax)
crosspairs(X,Y,rmax)
cobble.xy(x, y, f, fatal)
commasep(x)
contour.im(x, \dots)
countingweights(id, areas, check = TRUE)
damaged.ppm(object)
default.expand(object, m, epsilon)
default.clipwindow(object, epsilon)
default.n.tiling(X, nd, ntile, npix, verbose)
default.ntile(X)
densityhack(x, bw, adjust, kernel, window, width,
             give.Rkern, n, from, to, cut, na.rm, weights)
diagnose.ppm.engine(object, \dots, type="eem", typename, opt,
                         sigma=NULL, rbord = reach(object), compute.sd=TRUE,
                         compute.cts=TRUE, rv=NULL)
discretise(x, xrange, nx)
distpl(p, l)               
distppl(p, l)
distppll(p, l, mintype=0, method="Fortran", listit=FALSE)
divisors(n)
do.call.matched(fname, arglist, funargs, extrargs)
edge.Ripley(X, r, W, method)
edge.Trans(X, Y, W, exact, paired, trim)
ensure2vector(x)
equals.quad(Q)          
eratosthenes(nmax)
erode.mask(w,r)
even.breaks.owin(w)
exactdt(X, \dots)              
exactPdt(w)
fasp(fns, titles, formulae, which, dataname, title)
FormatFaspFormulae(f, argname)
fv(x, argu, ylab, valu, fmla, alim, labl, desc)
greatest.common.divisor(n,m)
getfields(X, L, fatal = TRUE)
getglmfit(object)
gridindex(x, y, xrange, yrange, nx, ny)            
handle.r.b.args(r = NULL, breaks = NULL, window, eps = NULL, rmaxdefault)
handle.rshift.args(W, \dots, radius, width, height, edge, clip, edgedefault)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
identical.formulae(x,y)
inside.xypolygon(pts, polly, test01 = TRUE, method="Fortran")     
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
is.cadlag(s)
is.data(Q)
is.fv(x)
is.marked(X, \dots)
is.marked.default(\dots)  
is.marked.ppp(X, na.action="warn", \dots)
is.poisson.ppm(x)
is.prime(n)
is.stationary.ppm(x)    
km.rs(o, cc, d, breaks)
Kborder.engine(X, rmax, nr, correction, weights)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks)
Kmulti.inhom(X, I, J, lambdaI, lambdaJ, \dots,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
         lambdaIJ=NULL,
         Iname = "points satisfying condition I",
         Jname = "points satisfying condition J")
killinteraction(model)
least.common.multiple(n,m)
\method{levels}{im}(x) <- value
lookup.im(Z, x, y, naok)
make.even.breaks(bmax, npos, bstep)
markspace.integral(X)
marks.quad(Q)           
matcolall(x)
matcolany(x)
matcolsum(x)            
matrixsample(mat, newdim, phase)
matrowall(x)
matrowany(x)
matrowsum(x)
maxflow(costm)
meanX.owin(w)            
meanY.owin(w)
mkcor(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
mpl.engine(Q, trend, interaction, \dots, covariates, correction,
	 rbord, use.gam, forcefit, callstring, precomputed, savecomputed)
mpl.get.covariates(covariates, locations, type)
mpl.prepare(Q, X, P, trend, interaction, covariates,
            want.trend, want.inter, correction, rbord, Pname,
            callstring, \dots, precomputed, savecomputed)
MultiPair.checkmatrix(mat, n, name)
nearest.pixel(x, y, im)
no.trend.ppm(x)
n.quad(Q)               
overlap.owin(A,B)
overlap.trapezium(xa, ya, xb, yb, verb = FALSE)
overlap.xypolygon(P, Q)
owinpolycheck(W, verbose=TRUE)
param.quad(Q)
pixellate(x, \dots, weights)
ploterodewin(W1, W2, col.edge, col.inside, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ)
plot.diagppm(x, \dots, which, plot.neg="image",
                         plot.smooth="imagecontour",
                         plot.sd=TRUE, spacing=0.1,
                         srange=NULL, monochrome=FALSE, main=NULL)
plot.minconfit(x, \dots)
plot.pppmatching(x, addmatch = NULL, main = NULL, \dots)
plot.profilepl(x, \dots, add=FALSE, main=NULL, tag=TRUE, coeff=NULL, xvariable=NULL)
plot.qqppm(x, \dots, limits=TRUE, monochrome=FALSE,
           limcol=if(monochrome) "black" else "red")
plot.quadratcount(x, \dots, add, entries, dx, dy)
plot.quadrattest(x, \dots)
polynom(x, \dots)
pppdist.prohorov(X, Y, precision = 7)
pppmatching(X, Y, am, ty = "generic")
primefactors(n, prmax)
print.diagppm(x, \dots)
print.envelope(x, \dots)
print.fasp(x, \dots)
print.fv(x, \dots)
print.interact(x, \dots)       
print.isf(x, \dots)
print.minconfit(x, \dots)
print.plotppm(x, \dots)
print.pppmatching(x, \dots)
print.profilepl(x, \dots)
print.qqppm(x, \dots)
print.rmhcontrol(x, \dots)
print.rmhmodel(x, \dots)
print.rmhstart(x, \dots)
print.rmhseed(x, \dots)
print.summary.owin(x, \dots)
print.summary.ppp(x, \dots, dp=3)
print.summary.psp(x, \dots)
quad(data, dummy, w, param)
quadrat.breaks(xr, yr, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
quadrat.count.engine(x, y, xbreaks, ybreaks, weights)
quadscheme.replicated(data, dummy, method = "grid", \dots)
quadscheme.spatial(data, dummy, method = "grid", \dots)
rasterx.im(x)
rastery.im(x)
rebound.ppp(x, \dots)
rebound.psp(x, \dots)
rebound.owin(w, rect)
resid4plot(RES, plot.neg="image", plot.smooth="imagecontour",
           spacing=0.1, srange=NULL,monochrome=FALSE, main=NULL, \dots)
resid1plot(RES, opt, plot.neg="image", plot.smooth="imagecontour",
              srange=NULL, monochrome=FALSE, main=NULL, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD, xlab,ylab, \dots)
resolve.defaults(\dots)
reverse.xypolygon(p)
rmax.rule(fun, W, lambda)
rotxy(X, angle = pi/2)
rmhResolveTypes(model, start, control)
rmhcontrol.rmhcontrol(control, \dots)
rmhcontrol.list(control, \dots)
rmhEngine(InfoList, \dots, verbose, reseed, kitchensink, preponly) 
rmhmodel.rmhmodel(model, \dots)
rmhmodel.list(model, \dots)
rmhseed(iseed)
rmhstart.rmhstart(start, \dots)
rmhstart.list(start, \dots)
rmhmodel.ppm(model, win, \dots, verbose, project, control)
rmpoint.I.allim(n, f, types)
rpoint.multi(n, f, fmax=NULL, marks = NULL, win = unit.square(), giveup = 1000, verbose = FALSE)
runifdisc(n, r = 1, x = 0, y = 0)
runifpoispp(lambda, win = owin(c(0, 1), c(0, 1)))
runifrect(n, win = owin(c(0, 1), c(0, 1)))
second.moment.calc(x, sigma, edge=TRUE, what="Kmeasure", debug=FALSE, \dots)
shiftxy(X, vec = c(0, 0))
spatstat.rawdata.location(\dots)             
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)             
stratrand(window, nx, ny, k = 1)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
summary.envelope(object,\dots)
summary.pppmatching(object, \dots)
sympoly(x, y, n)
termsinformula(x)
tilecentroids(W, nx, ny)        
trim.rectangle(W, xmargin, ymargin=xmargin)
update.interact(object, \dots)
validate.mask(w, fatal=TRUE)        
validate.quad(Q, fatal, repair, announce)
variablesinformula(x)
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
verify.xypolygon(p, fatal=TRUE)
whist(x,breaks,weights, trim, right)
w.quad(Q)               
x.quad(Q)
y.quad(Q)
xypolyselfint(p, eps, proper)
xypolygon2psp(p, w)
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
