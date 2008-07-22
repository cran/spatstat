subroutine badgey(u,v,ix,x,y,npts,par,period,cifval,ndisc,aux)
#
# Conditional intensity function for a Baddeley-Geyer process.
# Since this is an unmarked process the arguments
# 
#              mrk, marks, and ntypes
# 
# are dummy arguments.
#
implicit double precision (a-h,o-z)
integer aux(ndisc,1)
dimension x(1), y(1), marks(1), par(1), period(2)

# Hey --- you can ***do*** this, it seems! Dynamical
# memory allocation ... fortran 90 coming into play?
dimension w(ndisc), tee(ndisc)
logical newpt, per

beta  = par(1)
if(npts==0) {
	cifval = beta
	return
}

eps  = 2.22e-16 # Essentially .Machine$double.eps from Splus.
zero = 0.d0
one  = 1.d0
per  = period(1) > zero

if(ix > 0) {
	if(per) call dist2(u,v,x(ix),y(ix),period,d2)
	else d2 = (u-x(ix))**2 + (v-y(ix))**2
	newpt = d2 > eps
} else newpt = .true.

do i=1,ndisc {
	w(i) = zero
}

if(newpt) {
	do i=1,ndisc {
		tee(i) = zero
	}
	do j=1,npts {
		if(j == ix) next
		if(per) call dist2(u,v,x(j),y(j),period,d2)
		else d2 = (u-x(j))**2 + (v-y(j))**2
		do i=ndisc,1,-1 {
			r2    = par(3*i+1)
			s     = par(3*i+2)
			if(d2 < r2) {
				tee(i) = tee(i) + one
				a = aux(i,j)
				if(ix>0) { # Adjust
					if(per) call dist2(x(ix),y(ix),
                                                           x(j),y(j),period,dd2)
					else dd2 = (x(ix)-x(j))**2 +
                                                         (y(ix)-y(j))**2
					if(dd2 < r2) a = a - one
				}
				b = a + one
				if(a < s & s < b) {
					w(i) = w(i) + s - a
				}
				else if(s >= b) w(i) = w(i) + one
			} else break
		}
	}
} else {
	do i=1,ndisc {
		tee(i) = aux(i,ix)
	}
	do j=1,npts {
		if(j == ix) next
		if(per) call dist2(u,v,x(j),y(j),period,d2)
		else d2 = (u-x(j))**2 + (v-y(j))**2
		do i=ndisc,1,-1 {
			r2    = par(3*i+1)
			s     = par(3*i+2)
			if(d2 < r2) {
				a = aux(i,j) - 1
				b = a + one
				if(a < s & s < b) {
					w(i) = w(i) + s - a
				}
				else if(s >= b) w(i) = w(i) + one
			} else break
		}
	}
}
some = zero
do i=1,ndisc {
	s = par(3*i+2)
	w(i) = w(i) + min(tee(i),s)
	gamma = par(3*i)
	if(gamma < eps ) {
		if(w(i) > zero) {
			cifval = zero
			return
		} # else some = some + zero
	}
	else some = some + log(gamma)*w(i)
}
cifval = beta*exp(some)
return
end
