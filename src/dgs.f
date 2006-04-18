C Output from Public domain Ratfor, version 1.0
      subroutine dgs(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      dimension par(3), x(n), y(n), period(2)
      logical per
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      per = period(1) .gt. zero
      beta = par(1)
      rho = par(2)
      r2 = par(3)
      a = two*atan(one)/rho
      if(n.eq.0)then
      cifval = beta
      return
      endif
      soglum = zero
      do23002 j = 1,n 
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
23002 continue
23003 continue
      cifval = beta*exp(soglum)
      return
      end
