C Output from Public domain Ratfor, version 1.0
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
      do23000 j = 1,n 
      if(j .eq. ix)then
      sincr = zero
      else
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .lt. r2)then
      sincr = two*log(sin(a*sqrt(d2)))
      else
      sincr = zero
      endif
      endif
      soglum = soglum + sincr
23000 continue
23001 continue
      cifval = beta*exp(soglum)
      return
      end
