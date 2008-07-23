subroutine updaux(itype,x,y,u,v,npts,ix,par,period,ndisc,aux)
implicit double precision(a-h,o-z)
dimension x(1), y(1), par(2), period(2)
# Dimensioning par(2) is just a trick to keep g77 from bitching;
# par(1) is really adequate.
integer aux(ndisc,1)
logical per

ndisc = par(2)
per   = period(1) > 0.d0
nm1   = npts-1

if(itype == 1) { # Birth
	do k=1,ndisc {
		aux(k,npts) = 0
	}
	do j=1,nm1 {
		if(per) call dist2(u,v,x(j),y(j),period,d2)
		else d2 = (u-x(j))**2 + (v-y(j))**2
		do k=ndisc,1,-1 {
			r2 = par(3*k+1)
			if(d2 < r2) {
				aux(k,j) = aux(k,j)+1
				aux(k,npts) = aux(k,npts)+1
			} else break
		}
	}
	return
}

if(itype == 2) { # Death
	do j = 1,npts {
		if(j==ix) next
		if(per) call dist2(x(ix),y(ix),x(j),y(j),period,d2)
		else d2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
		do k=ndisc,1,-1 {
			r2 = par(3*k+1)
			if(d2 >= r2 & j < ix) next
			if(d2 < r2) {
				if(j < ix) aux(k,j) = aux(k,j) - 1
				else aux(k,j-1) = aux(k,j) - 1
			} else if(j>=ix)
					aux(k,j-1) = aux(k,j)
		}
	}
	do k = 1,ndisc {
		aux(k,npts) = 0
	}
	return
}

if(itype == 3) { # Shift
	do k = 1,ndisc {
		aux(k,ix) = 0
	}
	do j = 1,npts {
		if(j == ix) next
		if(per) call dist2(u,v,x(j),y(j),period,d2)
		else d2 = (u-x(j))**2 + (v-y(j))**2
		if(per) call dist2(x(ix),y(ix),x(j),y(j),period,d2i)
		else d2i = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
		do k = 1,ndisc {
			r2 = par(3*k+1)
			if(d2 >= r2 & d2i >= r2) next
			if(d2 < r2) {
				aux(k,ix) = aux(k,ix) + 1
				if(d2i >= r2) aux(k,j) = aux(k,j) + 1
			} else if(d2i < r2) aux(k,j) = aux(k,j) - 1
		}
	}
	return
}

call fexit("Argument itype to updaux must be 1, 2, or 3; bailing out.\n")
end
