C Output from Public domain Ratfor, version 1.0
      subroutine diggra(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      double precision kappa
      dimension par(7), x(n), y(n), period(2)
      logical per
      zero = 0.d0
      per = period(1) .gt. zero
      beta = par(1)
      kappa = par(2)
      delta = par(3)
      rho = par(4)
      d2 = par(5)
      r2 = par(6)
      a = par(7)
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
      call dist2(u,v,x(j),y(j),period,t2)
      else
      t2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(t2 .lt. d2)then
      cifval = zero
      return
      endif
      if(t2 .lt. r2)then
      sincr = log(sqrt(t2) - delta) - a
      else
      sincr = zero
      endif
      endif
      soglum = soglum + sincr
23002 continue
23003 continue
      cifval = beta*exp(kappa*soglum)
      return
      end
