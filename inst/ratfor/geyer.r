# $Id: geyer.r,v 1.4 2008/07/21 23:40:05 adrian Exp adrian $
subroutine geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
#
# Conditional intensity function for a Geyer process.
#

implicit double precision(a-h,o-z)
dimension par(5), x(1), y(1), period(2)
integer aux(1)
logical per, newpt

eps  = 2.22d-16 # Essentially .Machine$double.eps from Splus.
zero = 0.d0
one  = 1.d0
per  = period(1) > zero

beta  = par(1)
# ndisc = par(2) = 1
gamma = par(3)
r2    = par(4) # This is the ***squared*** radius; avoids needing square root.
s     = par(5)

if(npts==0) {
	cifval = beta
	return
}

# Check on type of event.  (Makes a difference to how the value
# of c1 = t((u,v),X) is calculated.):
#
# ix < 0          <--> proposing simple birth
# ix > 0 & !newpt <--> proposing death
# ix > 0 & newpt  <--> proposing shift from (x(ix), y(ix)) to (u, v).

if(ix > 0) {
	if(per) call dist2(u,v,x(ix),y(ix),period,d2)
	else d2 = (u-x(ix))**2 + (v-y(ix))**2
	newpt = d2 > eps
}
else newpt = .true.
w = zero
if(newpt) {
	tee = zero
	do j = 1,npts {
		if(j == ix) next
		if(per) call dist2(u,v,x(j),y(j),period,d2)
		else d2 = (u-x(j))**2 + (v-y(j))**2
		if(d2 < r2) {
			tee = tee + one
			a = aux(j)
			if(ix > 0) { # Adjust
				if(per) call dist2(x(ix),y(ix),x(j),y(j),period,dd2)
				else dd2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
				if(dd2 < r2) a = a - one
			}
			b = a + one
			if(a < s & s < b) {
				w = w + s - a
			}
			else if(s >= b) w = w + one
		}
	}
} else {
	tee = aux(ix)
	do j = 1,npts {
                if(j == ix) next
                if(per) call dist2(u,v,x(j),y(j),period,d2)
                else d2 = (u-x(j))**2 + (v-y(j))**2
                if(d2 < r2) {
                        a = aux(j) - 1
                        b = a + one
                        if(a < s & s < b) {
                                w = w + s - a
                        }
                        else if(s >= b) w = w + one
                }
	}

}
w = w + min(tee,s)
count = min(s,c1) + c2
if(gamma < eps ) {
	if(w > zero) cifval = zero
	else cifval = beta
}
else cifval = beta*exp(log(gamma)*w)

return
end
