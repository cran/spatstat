\name{spatstat-internal}
\title{Internal spatstat functions}
\alias{[.hyperframe}
\alias{[<-.hyperframe}
\alias{$.hyperframe}
\alias{$<-.hyperframe}
\alias{[.localpcfmatrix}
\alias{acedist.show}
\alias{acedist.noshow}
\alias{adjustthinrange}
\alias{affinexy}
\alias{affinexypolygon}
\alias{anycrossing.psp}
\alias{apply23sum}
\alias{area.xypolygon}       
\alias{areaGain.diri}       
\alias{areaGain.grid}       
\alias{areaLoss.diri}       
\alias{areaLoss.grid}       
\alias{AsymmDistance.psp}
\alias{as.breakpts}
\alias{as.data.frame.fv}
\alias{as.fv}
\alias{as.im.linim}
\alias{as.list.hyperframe}
\alias{as.listof}
\alias{as.units}
\alias{beachcolours}
\alias{bermantestEngine}
\alias{bdrylength.xypolygon} 
\alias{bdry.mask} 
\alias{blankcoefnames}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{bt.frame}
\alias{cartesian}
\alias{cat.factor}
\alias{cellmiddles}
\alias{censtimeCDFest}
\alias{checkfields}          
\alias{check.finite}
\alias{check.hist.lengths}
\alias{check.named.list}
\alias{check.named.thing}
\alias{check.named.vector}
\alias{check.nvector}
\alias{check.nmatrix}
\alias{check.1.real}
\alias{check.testfun}
\alias{clarkevansCalc}
\alias{clip.psp}
\alias{cliprect.psp}
\alias{clippoly.psp}
\alias{closepairs}
\alias{complaining}
\alias{compileK}
\alias{compilepcf}
\alias{crosspairs}
\alias{crosspaircounts}
\alias{crossfikselterms}
\alias{cobble.xy}
\alias{commasep}
\alias{compatible.units}
\alias{conform.imagelist}
\alias{countingweights}
\alias{damaged.ppm}
\alias{default.clipwindow}
\alias{default.n.tiling}
\alias{default.ntile}
\alias{densitypointsEngine}
\alias{dgsTerms}
\alias{diagnose.ppm.engine}
\alias{diggraterms}
\alias{digital.volume}
\alias{dilate.owin}
\alias{dim.fasp}               
\alias{dim.hyperframe}               
\alias{dim.im}               
\alias{dim.msr}
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{dimnames.msr}
\alias{distpl}               
\alias{distppl}
\alias{distppll}
\alias{distppllmin}
\alias{distributecbind}
\alias{dist2dpath}
\alias{divisors}
\alias{do.as.im}
\alias{do.call.matched}
\alias{do.iplot}
\alias{do.istat}
\alias{edge.Ripley}
\alias{edge.Trans}
\alias{edge.Trans.modif}
\alias{emptywindow}
\alias{ensure2vector}
\alias{envelopeEngine}
\alias{equalpairs}          
\alias{equalpairs.quad}          
\alias{equals.quad}          
\alias{equalsfun.quad}          
\alias{eratosthenes}
\alias{erodemask}
\alias{erode.owin}
\alias{evalCovar}
\alias{evalCovariate}
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{even.breaks.owin}
\alias{evenly.spaced}
\alias{exactdt}              
\alias{exactPdt}
\alias{explain.ifnot}
\alias{extractAIC.slrm}
\alias{f3engine}
\alias{f3Cengine}
\alias{fasp}
\alias{findCovariate}
\alias{findcbind}
\alias{findmarktype}
\alias{firstfactor}
\alias{fii}
\alias{fikselterms}
\alias{fillNA}
\alias{flipxypolygon}
\alias{forbidNA}
\alias{FormatFaspFormulae}
\alias{formula.slrm}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{fvnames}
\alias{fvnames<-}
\alias{g3engine}
\alias{g3Cengine}
\alias{gpcmethod}
\alias{greatest.common.divisor}
\alias{getdataname}
\alias{getfields}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{geyercounts}
\alias{GLMpredict}
\alias{gridindex}            
\alias{grid1index}
\alias{grow.rectangle}
\alias{handle.r.b.args}
\alias{handle.rshift.args}
\alias{ho.engine}
\alias{identical.formulae}
\alias{idorempty}
\alias{implemented.for.K}
\alias{inpoint}
\alias{inside.triangle}     
\alias{inside.xypolygon}     
\alias{intersect.ranges}     
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}     
\alias{is.cadlag}
\alias{is.data}
\alias{is.fv}
\alias{is.hole.xypolygon}
\alias{is.hyperframe}
\alias{is.infline}
\alias{is.interact}
\alias{is.marked.default}    
\alias{is.marked.psp}    
\alias{is.marked.quad}    
\alias{is.multitype.quad}    
\alias{is.multitype.default}    
\alias{is.parseable}
\alias{is.pp3}
\alias{is.ppx}
\alias{is.prime}
\alias{is.psp}
\alias{is.slrm}
\alias{is.tess}
\alias{k3engine}
\alias{Kborder.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{Kmulti.inhom}               
\alias{Kpcf.kppm}               
\alias{killinteraction}
\alias{km.rs.opt}
\alias{least.common.multiple}
\alias{levels.im}
\alias{levels<-.im}
\alias{licence.polygons}
\alias{license.polygons}
\alias{linequad}
\alias{linearKengine}
\alias{linearpcfengine}
\alias{listof}
\alias{localKengine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{paste.expr}
\alias{print.localpcfmatrix}
\alias{plot.localpcfmatrix}
\alias{lookup.im}
\alias{majorminorversion}
\alias{make.even.breaks}
\alias{make.parseable}
\alias{makeunits}
\alias{markappend}
\alias{markcbind}
\alias{markformat}
\alias{markformat.ppp}
\alias{markformat.psp}
\alias{markformat.default}
\alias{mark.scale.default}
\alias{markspace.integral}
\alias{marks.default}           
\alias{marks.quad}           
\alias{\%mapp\%}
\alias{markappendop}
\alias{marksubset}
\alias{markreplicateop}
\alias{\%mrep\%}
\alias{marksubsetop}
\alias{\%msub\%}
\alias{mask2df}
\alias{matcolall}
\alias{matcolany}
\alias{matcolsum}
\alias{matrixsample}
\alias{matrowall}
\alias{matrowany}
\alias{matrowsum}
\alias{maxflow}
\alias{meanlistfv}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{sensiblevarname}
\alias{sewpcf}
\alias{sewsmod}
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{model.matrix.slrm}
\alias{model.se.image}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{MultiPair.checkmatrix}
\alias{multiply.only.finite.entries}
\alias{na.handle.im}
\alias{names.hyperframe}
\alias{names<-.hyperframe}
\alias{nearest.pixel}
\alias{nearest.valid.pixel}
\alias{newstyle.coeff.handling}
\alias{nncleanEngine}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{nobjects}
\alias{nobjects.ppp}
\alias{nobjects.ppx}
\alias{nobjects.psp}
\alias{n.quad}
\alias{numalign}
\alias{nzpaste}
\alias{onecolumn}
\alias{ordinal}
\alias{outdated.interact}
\alias{overlap.owin}
\alias{overlap.trapezium}
\alias{overlap.xypolygon}    
\alias{oversize.quad}    
\alias{owinpolycheck}
\alias{owinpoly2mask}
\alias{pairs.listof}
\alias{param.quad}
\alias{paren}
\alias{partialModelMatrix}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{pickoption}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.barplotdata}
\alias{plot.bw.relrisk}
\alias{plot.fii}
\alias{plot.minconfit}
\alias{plot.plotpairsim}
\alias{plot.pppmatching}
\alias{plot.profilepl}
\alias{plot.qqppm}
\alias{plot.quadratcount}
\alias{plot.quadrattest}
\alias{polynom}
\alias{ppllengine}
\alias{ppmCovariates}
\alias{pppdist.mat}
\alias{pppdist.prohorov}
\alias{prange}
\alias{prefixfv}
\alias{primefactors}       
\alias{print.bt.frame}
\alias{print.bw.relrisk}
\alias{print.colourmap}
\alias{print.diagppm}
\alias{print.envelope}
\alias{print.fasp}       
\alias{print.fii}       
\alias{print.fv}       
\alias{print.hyperframe}
\alias{print.interact}       
\alias{print.isf}
\alias{print.linim}
\alias{print.lut}
\alias{print.minconfit}
\alias{print.msr}
\alias{print.plotpairsim}
\alias{print.plotppm}
\alias{print.pppmatching}
\alias{print.profilepl}
\alias{print.quadrattest}
\alias{print.qqppm}
\alias{print.rmhcontrol}
\alias{print.rmhmodel}
\alias{print.rmhstart}
\alias{print.slrm}
\alias{print.splitppp}
\alias{print.summary.fii}
\alias{print.summary.hyperframe}
\alias{print.summary.listof}
\alias{print.summary.lut}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{print.summary.psp}
\alias{print.summary.splitppp}
\alias{print.summary.units}
\alias{print.tess}
\alias{print.units}
\alias{prolongseq}
\alias{quad}
\alias{rectquadrat.breaks}
\alias{rectquadrat.countEngine}
\alias{rhs.of.formula}
\alias{quadrat.testEngine}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{pointgrid}
\alias{rasterx.im}
\alias{rastery.im}
\alias{reach.fii}
\alias{rebadge.fv}
\alias{rebound}
\alias{rebound.im}
\alias{rebound.ppp}
\alias{rebound.psp}
\alias{rebound.owin}
\alias{reconcile.fv}
\alias{repair.old.factor.image}
\alias{reincarnate.interact}
\alias{rescale.units}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{resolve.defaults}
\alias{resolve.2D.kernel}
\alias{restrict.mask}
\alias{reverse.xypolygon}    
\alias{revcumsum}
\alias{rmax.rule}
\alias{rotxy}
\alias{rotxypolygon}
\alias{row.names.hyperframe}
\alias{row.names<-.hyperframe}
\alias{runifpoispp}          
\alias{runifpoisppOnLines}          
\alias{runifrect}
\alias{rmhResolveExpansion}
\alias{rmhResolveTypes}
\alias{rmhcontrol.rmhcontrol}
\alias{rmhcontrol.list}
\alias{rmhEngine}
\alias{rmhmodel.rmhmodel}
\alias{rmhstart.rmhstart}
\alias{rmhstart.list}
\alias{rmpoint.I.allim}
\alias{rpoint.multi}
\alias{safelookup}
\alias{second.moment.calc}
\alias{second.moment.engine}
\alias{shiftxy}              
\alias{shiftxypolygon}              
\alias{simplify.xypolygon}              
\alias{simulrecipe}              
\alias{singlestring}              
\alias{smoothpointsEngine}              
\alias{sort.im}
\alias{spatstat.rawdata.location}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{sphere.volume}
\alias{stratrand}
\alias{strausscounts}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summary.envelope}
\alias{summary.fii}
\alias{summary.hyperframe}
\alias{summary.lut}
\alias{summary.profilepl}
\alias{summary.pppmatching}
\alias{summary.ppx}
\alias{summary.units}
\alias{superimposeMarks}
\alias{sympoly}
\alias{terms.slrm}
\alias{termsinformula}
\alias{test.crossing.psp}
\alias{test.selfcrossing.psp}
\alias{tilecentroids}        
\alias{tilenames}        
\alias{trim.mask}        
\alias{tweak.fv.entry}        
\alias{unitname.default}
\alias{unitname<-.default}
\alias{update.interact}
\alias{update.slrm}
\alias{validradius}        
\alias{valid.ppm}        
\alias{validate.mask}        
\alias{validate.quad}        
\alias{validposint}        
\alias{variablesinformula}
\alias{vcov.slrm}
\alias{verbalogic}
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{versionstring.spatstat}
\alias{verifyclass}
\alias{verify.xypolygon}
\alias{warn.ignored.args}
\alias{with.msr}
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}
\alias{xtfrm.im}
\alias{xypolygon2psp}
\alias{xypolyselfint}
\description{
  Internal spatstat functions.
}
\usage{
\method{[}{hyperframe}(x, i, j, drop=FALSE, ...)
\method{[}{hyperframe}(x, i, j) <- value
\method{$}{hyperframe}(x, name)
\method{$}{hyperframe}(x, i) <- value
acedist.show(X, Y, n, d, timelag = 0)
acedist.noshow(X, Y, n, d)
adjustthinrange(ur,vstep,vr)
affinexy(X, mat = diag(c(1, 1)), vec = c(0, 0), invert=FALSE)
affinexypolygon(p, mat, vec, detmat)
anycrossing.psp(A,B)
apply23sum(x)
area.xypolygon(polly)
areaGain.diri(u, X, r, ..., W=as.owin(X))
areaGain.grid(u, X, r, ..., W=NULL, ngrid=spatstat.options("ngrid.disc"))
areaLoss.diri(X, r, ..., W=as.owin(X), subset=NULL)
areaLoss.grid(X, r, ..., W=as.owin(X), subset=NULL, ngrid=spatstat.options("ngrid.disc"))
AsymmDistance.psp(X, Y, metric="Hausdorff", method="Fortran")
as.breakpts(\dots)
\method{as.data.frame}{fv}(x, ...)
as.fv(x)
\method{as.list}{hyperframe}(x, ...)
\method{as.im}{linim}(X, \dots)
as.listof(x)
as.units(s)
beachcolours(heightrange, sealevel, monochrome, ncolours)
bermantestEngine(model, covariate, which, alternative, ..., modelname, covname, dataname)
bdrylength.xypolygon(polly)
bdry.mask(W)
blankcoefnames(x)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
bt.frame(Q, trend=~1, interaction=NULL, ..., covariates=NULL,
         correction="border", rbord=0, use.gam=FALSE, allcovar=FALSE)
cartesian(pp, markset, fac = TRUE)
cat.factor(\dots, recursive=FALSE)
cellmiddles(W, nx, ny, npix, gi)
censtimeCDFest(o, cc, d, breaks, \dots, KM, RS, HAN, RAW, han.denom, pmax)
checkfields(X,L)          
check.finite(x, context, xname)
check.hist.lengths(hist,breaks)
check.named.list(x, nam, context, namopt)
check.named.vector(x, nam, context, namopt)
check.named.thing(x, nam, namopt, xtitle, valid, type, context, fatal)
check.nvector(v, npoints, fatal=TRUE, things="data points", naok=FALSE)
check.nmatrix(m, npoints, fatal=TRUE, things="data points", naok=FALSE, squarematrix=TRUE)
check.1.real(x, context, fatal)
check.testfun(f, f1, X)
clarkevansCalc(X, correction, clipregion, working)
clip.psp(x, window, check=TRUE)
cliprect.psp(x, window)
clippoly.psp(s, window)
closepairs(X,rmax)
complaining(whinge, fatal, value)
compileK(D, r, weights, denom, check)
compilepcf(D, r, weights, denom, check, endcorrect, \dots)
crosspairs(X,Y,rmax)
crosspaircounts(X,Y,r)
crossfikselterms(X,Y,r,kappa)
cobble.xy(x, y, f, fatal)
commasep(x)
compatible.units(x, y, coerce=TRUE)
conform.imagelist(X, Zlist)
countingweights(id, areas, check = TRUE)
damaged.ppm(object)
default.clipwindow(object, epsilon)
default.n.tiling(X, nd, ntile, npix, verbose)
default.ntile(X)
densitypointsEngine(x, sigma, \dots, weights, edge, varcov,leaveoneout,diggle)
diagnose.ppm.engine(object, \dots, type="eem", typename, opt,
                         sigma=NULL, rbord = reach(object), compute.sd=TRUE,
                         compute.cts=TRUE, rv=NULL, oldstyle=FALSE)
dgsTerms(X, Y, idX, idY, rho)
diggraterms(X, Y, idX, idY, delta, rho)
digital.volume(range, nval, vside)
dilate.owin(\dots)
\method{dim}{fasp}(x)
\method{dim}{hyperframe}(x)
\method{dim}{im}(x)
\method{dim}{msr}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
\method{dimnames}{msr}(x)
distpl(p, l)               
distppl(p, l)
distppll(p, l, mintype=0, method="Fortran", listit=FALSE)
distppllmin(p, l)
distributecbind(x)
dist2dpath(dist, method="C")
divisors(n)
do.as.im(x, action, \dots, W, eps, dimyx, xy, na.replace)
do.call.matched(fname, arglist, funargs, extrargs)
do.iplot(panel)
do.istat(panel)
edge.Ripley(X, r, W, method)
edge.Trans(X, Y, W, exact, paired, trim)
edge.Trans.modif(X, Y, WX, WY, exact, paired, trim)
emptywindow(w)
ensure2vector(x)
envelopeEngine(X, fun, simul,
           nsim=99, nrank=1, \dots, 
           verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL, maxnerr=nsim, internal=NULL, cl=NULL,
           envir.user=envir.user)
equalpairs(U, X, marked=FALSE)
equalpairs.quad(Q)
equals.quad(Q)          
equalsfun.quad(Q)          
eratosthenes(nmax)
erodemask(w,r,strict)
erode.owin(\dots)
evalCovar(model, covariate, \dots,
          dimyx, eps, jitter, modelname, covname, dataname)
evalCovariate(covariate, locations)
evalInteraction(X,P,E,interaction,correction,\dots,precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,\dots,precomputed,savecomputed)
even.breaks.owin(w)
evenly.spaced(x, tol)
exactdt(X, \dots)              
exactPdt(w)
explain.ifnot(expr, context)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, ...)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
fasp(fns, which, formulae, dataname, title, rowNames, colNames)
findCovariate(covname, scope, scopename=NULL)
findcbind(root, depth, maxdepth)
findmarktype(x)
firstfactor(x)
fii(interaction, coefs, Vnames, IsOffset)
fikselterms(U,X,r,kappa,EqualPairs)
fillNA(x, value)
flipxypolygon(p)
forbidNA(x, context, xname)
FormatFaspFormulae(f, argname)
\method{formula}{slrm}(x, ...)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
fvnames(X, a)
fvnames(X, a) <- value
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
gpcmethod(fname, signature)
greatest.common.divisor(n,m)
getdataname(defaultvalue, \dots, dataname)
getfields(X, L, fatal = TRUE)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
GLMpredict(fit, data, coefs, changecoef)
gridindex(x, y, xrange, yrange, nx, ny)            
grid1index(x, xrange, nx)
grow.rectangle(W, xmargin=0, ymargin=xmargin)
handle.r.b.args(r = NULL, breaks = NULL, window, eps = NULL, rmaxdefault)
handle.rshift.args(W, \dots, radius, width, height, edge, clip, edgedefault)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
identical.formulae(x,y)
idorempty(w, r, caller)
implemented.for.K(correction, windowtype, explicit)
inpoint(W)
inside.triangle(x, y, xx, yy)
inside.xypolygon(pts, polly, test01 = TRUE, method="Fortran")
intersect.ranges(a,b,fatal)
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
is.cadlag(s)
is.data(Q)
is.fv(x)
is.hole.xypolygon(polly)
is.hyperframe(x)
is.infline(x)
is.interact(x)
\method{is.marked}{default}(\dots)  
\method{is.marked}{psp}(X, \dots)
\method{is.marked}{quad}(X, na.action="warn", \dots)
\method{is.multitype}{default}(\dots)  
\method{is.multitype}{quad}(X, na.action="warn", \dots)
is.parseable(x)
is.pp3(x)
is.ppx(x)
is.prime(n)
is.psp(x)
is.slrm(x)
is.tess(x)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks)
Kmulti.inhom(X, I, J, lambdaI=NULL, lambdaJ=NULL, \dots,
         r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
         lambdaIJ=NULL, sigma=NULL, varcov=NULL, 
         Iname = "points satisfying condition I",
         Jname = "points satisfying condition J")
Kpcf.kppm(model, what)
killinteraction(model)
km.rs.opt(o, cc, d, breaks, KM, RS)
least.common.multiple(n,m)
\method{levels}{im}(x)
\method{levels}{im}(x) <- value
licence.polygons()
license.polygons()
linequad(X, Y, ..., eps, nd)
linearKengine(X, \dots, r, reweight, denom, correction, showworking)
linearpcfengine(X, \dots, r, reweight, denom, correction)
listof(\dots)
localKengine(X, \dots, wantL, lambda, correction, verbose, rvalue)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
paste.expr(x)
\method{print}{localpcfmatrix}(x, \dots)
\method{plot}{localpcfmatrix}(x, \dots)
\method{[}{localpcfmatrix}(x, i, ...)
lookup.im(Z, x, y, naok, strict)
majorminorversion(v)
make.even.breaks(bmax, npos, bstep)
make.parseable(x)
makeunits(sing, plur, mul)
markappend(\dots)
markcbind(\dots)
markformat(x)
\method{markformat}{ppp}(x)
\method{markformat}{psp}(x)
\method{markformat}{default}(x)
mark.scale.default(marx, w, markscale=NULL, maxsize=NULL)
markspace.integral(X)
\method{marks}{default}(x, \dots)
\method{marks}{quad}(x, dfok=FALSE, \dots)
markappendop(x, y)
marksubset(x, index, format)
marksubsetop(x, i)
markreplicateop(x, n)
mask2df(w)
matcolall(x)
matcolany(x)
matcolsum(x)            
matrixsample(mat, newdim, phase)
matrowall(x)
matrowany(x)
matrowsum(x)
maxflow(costm)
meanlistfv(z)
meanX.owin(w)            
meanY.owin(w)
sensiblevarname(guess, fallback, maxlen)
sewpcf(d, w, denargs, lambda2area)
sewsmod(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
\method{model.matrix}{slrm}(object, ..., keepNA = TRUE)
model.se.image(fit, W, \dots, what)
mpl.engine(Q, trend, interaction, \dots, covariates, covfunargs, correction,
	 rbord, use.gam, gcontrol, famille,
         forcefit, nd, allcovar, callstring, precomputed, savecomputed, preponly)
mpl.get.covariates(covariates, locations, type, covfunargs)
mpl.prepare(Q, X, P, trend, interaction, covariates, 
            want.trend, want.inter, correction, rbord, Pname,
            callstring, \dots, covfunargs, allcovar, precomputed, savecomputed,
            vnamebase, vnameprefix)
MultiPair.checkmatrix(mat, n, matname, naok, zerook)
multiply.only.finite.entries(x, a)
na.handle.im(X, na.replace)
\method{names}{hyperframe}(x)
\method{names}{hyperframe}(x) <- value
nearest.pixel(x, y, im)
nearest.valid.pixel(x, y, im)
newstyle.coeff.handling(object)
nncleanEngine(kthNND, k, d, \dots, tol, plothist, verbose, maxit)
nndcumfun(X, \dots, r)
no.trend.ppm(x)
nobjects(x)
\method{nobjects}{ppp}(x)
\method{nobjects}{ppx}(x)
\method{nobjects}{psp}(x)
n.quad(Q)
numalign(i, nmax, zero)
nzpaste(\dots, sep, collapse)
onecolumn(m)
ordinal(k)
outdated.interact(object)
overlap.owin(A,B)
overlap.trapezium(xa, ya, xb, yb, verb = FALSE)
overlap.xypolygon(P, Q)
oversize.quad(Q, \dots, nU, nX)
owinpolycheck(W, verbose=TRUE)
owinpoly2mask(w, rasta, check=TRUE)
\method{pairs}{listof}(\dots, plot=TRUE)
param.quad(Q)
paren(x, type)
partialModelMatrix(X,D,model,callstring,\dots)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI = NULL, lambdaJ = NULL, ...,
               r = NULL, breaks = NULL,
               kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
               correction = c("translate", "Ripley"),
               sigma = NULL, varcov = NULL,
               Iname = "points satisfying condition I",
               Jname = "points satisfying condition J")
pickoption(what="option", key, keymap, \dots,
           exact=FALSE, list.on.err=TRUE, die=TRUE, multi=FALSE)
ploterodewin(W1, W2, col.edge, col.inside, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ)
\method{plot}{barplotdata}(x, \dots)
\method{plot}{bw.relrisk}(x, \dots)
\method{plot}{fii}(x, \dots)
\method{plot}{minconfit}(x, \dots)
\method{plot}{pppmatching}(x, addmatch = NULL, main = NULL, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{profilepl}(x, \dots, add=FALSE, main=NULL, tag=TRUE, coeff=NULL, xvariable=NULL)
\method{plot}{qqppm}(x, \dots, limits=TRUE, monochrome=FALSE,
           limcol=if(monochrome) "black" else "red")
\method{plot}{quadratcount}(x, \dots, add, entries, dx, dy, show.tiles)
\method{plot}{quadrattest}(x, \dots)
polynom(x, \dots)
ppllengine(X, Y, action="project", check=FALSE)
ppmCovariates(model)
pppdist.mat(X, Y, cutoff = 1, q = 1, matching = TRUE, precision = 9, approximation = 10)
pppdist.prohorov(X, Y, n, dfix, type, cutoff = 1, matching = TRUE, ccode = TRUE, precision = 9, approximation = 10) 
prange(x)
prefixfv(x, prefix, prephrase)
primefactors(n, prmax)
\method{print}{bt.frame}(x, \dots)
\method{print}{bw.relrisk}(x, \dots)
\method{print}{colourmap}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{fii}(x, \dots)
\method{print}{fv}(x, \dots)
\method{print}{hyperframe}(x, ...)
\method{print}{interact}(x, \dots, family=TRUE, brief=FALSE)       
\method{print}{isf}(x, \dots)
\method{print}{linim}(x, \dots)
\method{print}{lut}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{plotppm}(x, \dots)
\method{print}{plotpairsim}(x, \dots)
\method{print}{pppmatching}(x, \dots)
\method{print}{profilepl}(x, \dots)
\method{print}{quadrattest}(x, \dots)
\method{print}{qqppm}(x, \dots)
\method{print}{rmhcontrol}(x, \dots)
\method{print}{rmhmodel}(x, \dots)
\method{print}{rmhstart}(x, \dots)
\method{print}{slrm}(x, ...)
\method{print}{splitppp}(x, \dots)
\method{print}{summary.fii}(x, \dots)
\method{print}{summary.hyperframe}(x, ...)
\method{print}{summary.listof}(x, \dots)
\method{print}{summary.owin}(x, \dots)
\method{print}{summary.ppp}(x, \dots, dp=3)
\method{print}{summary.psp}(x, \dots)
\method{print}{summary.splitppp}(x, \dots)
\method{print}{summary.units}(x, \dots)
\method{print}{tess}(x, \dots, brief=FALSE)
\method{print}{units}(x, \dots)
prolongseq(x, newrange)
quad(data, dummy, w, param)
rectquadrat.breaks(xr, yr, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
rectquadrat.countEngine(x, y, xbreaks, ybreaks, weights)
rhs.of.formula(x)
quadrat.testEngine(X, nx, ny, \dots, Xcount, xbreaks, ybreaks, tess=NULL, fit=NULL, Xname=NULL, fitname=NULL)
quadscheme.replicated(data, dummy, method = "grid", \dots)
quadscheme.spatial(data, dummy, method = "grid", \dots)
pointgrid(W, ngrid)
rasterx.im(x)
rastery.im(x)
\method{reach}{fii}(x, \dots, epsilon)
rebadge.fv(x, new.ylab, new.fname, tags, new.desc, new.labl, new.yexp,
new.dotnames, new.preferred, new.formula)
rebound(x, rect)
\method{rebound}{im}(x, rect)
\method{rebound}{ppp}(x, rect)
\method{rebound}{psp}(x, rect)
\method{rebound}{owin}(x, rect)
reconcile.fv(\dots)
repair.old.factor.image(x)
reincarnate.interact(object)
\method{rescale}{units}(X,s)
resid4plot(RES, plot.neg="image", plot.smooth="imagecontour",
           spacing=0.1, srange=NULL,monochrome=FALSE, main=NULL, \dots)
resid1plot(RES, opt, plot.neg="image", plot.smooth="imagecontour",
              srange=NULL, monochrome=FALSE, main=NULL, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD, xlab,ylab, \dots)
resolve.defaults(\dots, .StripNull=FALSE)
resolve.2D.kernel(\dots, sigma, varcov, x, mindist, adjust, useN)
restrict.mask(M, W)
reverse.xypolygon(p, adjust=FALSE)
revcumsum(x)
rmax.rule(fun, W, lambda)
rotxy(X, angle = pi/2)
rotxypolygon(p, angle = pi/2)
rmhResolveExpansion(win, control, imagelist, itype)
rmhResolveTypes(model, start, control)
\method{rmhcontrol}{rmhcontrol}(\dots)
\method{rmhcontrol}{list}(\dots)
rmhEngine(InfoList, \dots, verbose, track, kitchensink, preponly) 
\method{rmhmodel}{rmhmodel}(model, \dots)
\method{rmhstart}{rmhstart}(start, \dots)
\method{rmhstart}{list}(start, \dots)
rmpoint.I.allim(n, f, types)
\method{row.names}{hyperframe}(x)
\method{row.names}{hyperframe}(x) <- value
rpoint.multi(n, f, fmax=NULL, marks = NULL, win = unit.square(), giveup = 1000, verbose = FALSE)
runifpoispp(lambda, win = owin(c(0, 1), c(0, 1)))
runifpoisppOnLines(lambda, L)
runifrect(n, win = owin(c(0, 1), c(0, 1)))
safelookup(Z, X, factor, warn)
second.moment.calc(x, sigma=NULL, edge=TRUE, what="Kmeasure", debug=FALSE,
\dots, varcov=NULL, expand=FALSE)
second.moment.engine(x, sigma=NULL, edge=TRUE, what="Kmeasure", debug=FALSE,
\dots, obswin=as.owin(x), varcov=NULL, npts=NULL)
shiftxy(X, vec = c(0, 0))
shiftxypolygon(p, vec = c(0, 0))
simplify.xypolygon(p, dmin)
simulrecipe(type, expr, envir, csr)
singlestring(s, coll)
smoothpointsEngine(x, values, sigma, \dots, weights, varcov, leaveoneout)
\method{sort}{im}(x, \dots)
spatstat.rawdata.location(\dots)
spatialCDFframe(model, covariate, \dots)
spatialCDFtest(model, covariate, test, \dots, dimyx, eps, jitter, modelname, covname, dataname)
sphere.volume(range, nval = 10)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
stratrand(window, nx, ny, k = 1)
strausscounts(U,X,r,EqualPairs)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
\method{summary}{envelope}(object,\dots)
\method{summary}{fii}(object,\dots)
\method{summary}{hyperframe}(object, ..., brief=FALSE)
\method{summary}{profilepl}(object, \dots)
\method{summary}{pppmatching}(object, \dots)
\method{summary}{ppx}(object, \dots)
\method{summary}{units}(object, \dots)
superimposeMarks(arglist, nobj)
sympoly(x, y, n)
\method{terms}{slrm}(x, ...)
termsinformula(x)
test.crossing.psp(A,B)
test.selfcrossing.psp(A)
tilecentroids(W, nx, ny)
tilenames(x)
trim.mask(M, R, tolerant)
tweak.fv.entry(x, current.tag, new.labl=NULL, new.desc=NULL, new.tag=NULL)
\method{unitname}{default}(x)
\method{unitname}{default}(x) <- value
\method{update}{interact}(object, \dots)
\method{update}{slrm}(object, ..., evaluate = TRUE, env = parent.frame())
validradius(r, caller)
valid.ppm(object, na.value)
validate.mask(w, fatal=TRUE)        
validate.quad(Q, fatal, repair, announce)
validposint(n, caller, fatal)
variablesinformula(x)
\method{vcov}{slrm}(object, ...)
verbalogic(x, op)
versionstring.interact(object)
versionstring.ppm(object)
versionstring.spatstat()
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
verify.xypolygon(p, fatal=TRUE)
warn.ignored.args(\dots, context)
\method{with}{msr}(data, expr, \dots)
w.quad(Q)               
x.quad(Q)
y.quad(Q)
\method{xtfrm}{im}(x)
xypolyselfint(p, eps, proper, yesorno, checkinternal)
xypolygon2psp(p, w, check)
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}