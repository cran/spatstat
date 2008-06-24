require(spatstat)
source("marcelinodata.R")

I.lam <- predict (ppm(I.ppp, ~polynom(x,y,2)), type="trend")
J.lam <- predict (ppm(J.ppp, ~polynom(x,y,2)), type="trend")


Kinhom(I.ppp, lambda=I.lam, correction="iso")
Kinhom(I.ppp, lambda=I.lam, correction="border")

IJ.ppp=superimpose(a=I.ppp, b=J.ppp)

Kcross.inhom(IJ.ppp, i="a", j="b", I.lam, J.lam)
Kcross.inhom(IJ.ppp, i="a", j="b", I.lam, J.lam, correction = "iso")
Kcross.inhom(IJ.ppp, i="a", j="b", I.lam, J.lam, correction="border")


