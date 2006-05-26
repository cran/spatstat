C Output from Public domain Ratfor, version 1.0
      subroutine strauss(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      dimension par(3), x(1), y(1), period(2)
      logical per
      eps = 2.22d-16
      zero = 0.d0
      per = period(1) .gt. zero
      beta = par(1)
      gamma = par(2)
      r = par(3)
      if(n.eq.0)then
      cifval = beta
      return
      endif
      kount = 0
      if(per)then
      ixm1 = ix - 1
      ixp1 = max(1, ix + 1)
      if(ixm1 .gt. 0)then
      do23006 j = 1,ixm1 
      call dist2(u,v,x(j),y(j),period,d2)
      if(d2 .lt. r)then
      kount = kount+1
      endif
23006 continue
23007 continue
      endif
      if(ixp1 .le. n)then
      do23012 j = ixp1,n 
      call dist2(u,v,x(j),y(j),period,d2)
      if(d2 .lt. r)then
      kount = kount+1
      endif
23012 continue
23013 continue
      endif
      else
      ixm1 = ix - 1
      ixp1 = max(1, ix + 1)
      if(ixm1 .gt. 0)then
      do23018 j = 1,ixm1 
      a = r - (u - x(j))**2
      if(a .gt. 0)then
      a = a - (v-y(j))**2
      if(a .gt. 0)then
      kount = kount+1
      endif
      endif
23018 continue
23019 continue
      endif
      if(ixp1 .le. n)then
      do23026 j = ixp1,n 
      a = r - (u - x(j))**2
      if(a .gt. 0)then
      a = a - (v-y(j))**2
      if(a .gt. 0)then
      kount = kount+1
      endif
      endif
23026 continue
23027 continue
      endif
      endif
      if(gamma .lt. eps )then
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
