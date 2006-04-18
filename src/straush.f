C Output from Public domain Ratfor, version 1.0
      subroutine straush(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      dimension par(4), x(n), y(n), period(2)
      logical per
      eps = 2.22d-16
      zero = 0.d0
      per = period(1) .gt. zero
      beta = par(1)
      gamma = par(2)
      ri = par(3)
      rhc = par(4)
      if(n.eq.0)then
      cifval = beta
      return
      endif
      kount = 0
      do23002 j = 1,n 
      if(j .eq. ix)then
      continue
      else
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .lt. rhc)then
      cifval = zero
      return
      endif
      if(d2 .lt. ri)then
      kount = kount+1
      endif
      endif
23002 continue
23003 continue
      if(gamma .lt. eps)then
      if(kount .gt. 0)then
      cifval = zero
      else
      cifval = beta
      endif
      else
      cifval = beta*exp(log(gamma)*dble(kount))
      endif
      return
      end
