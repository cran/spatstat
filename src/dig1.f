      subroutine dig1(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      dimension par(2), x(n), y(n), period(2)
      logical per
      zero = 0.d0
      two = 2.d0
      per = period(1) .gt. 0.d0
      beta = par(1)
      rho = par(2)
      r2 = rho**2
      a = two*atan(1.d0)/rho
      soglum = zero
      do 23000 j = 1,n 
      if(.not.(j .eq. ix))goto 23002
      sincr = zero
      goto 23003
23002 continue
      if(.not.(per))goto 23004
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23005
23004 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23005 continue
      if(.not.(d2 .lt. r2))goto 23006
      sincr = two*log(sin(a*sqrt(d2)))
      goto 23007
23006 continue
      sincr = zero
23007 continue
23003 continue
      soglum = soglum + sincr
23000 continue
      cifval = beta*exp(soglum)
      return
      end
