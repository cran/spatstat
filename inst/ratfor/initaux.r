subroutine initaux(nmbr,par,period,x,y,npts,ndisc,aux)
implicit double precision(a-h,o-z)
dimension x(1), y(1), par(1), period(2)
integer aux(ndisc,1)
logical per

if(nmbr != 8 & nmbr != 11) return
per   = period(1) > 0.d0

do i=1,npts {
	do k=1,ndisc {
		aux(k,i) = 0
	}
	do j=1,npts {
		if(j==i) next
		if(per) call dist2(x(i),y(i),x(j),y(j),period,d2)
		else d2 = (x(i)-x(j))**2 + (y(i)-y(j))**2
		do k = ndisc,1,-1 {
			r2 = par(3*k+1)
			if(d2 < r2) {
				aux(k,i) = aux(k,i) + 1
			}
			else break
		}
	}
}

return
end
