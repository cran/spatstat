C Output from Public domain Ratfor, version 1.0
      subroutine geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
      implicit double precision(a-h,o-z)
      dimension par(4), x(1), y(1), period(2)
      integer aux(1)
      logical per, newpt
      eps = 2.22d-16
      zero = 0.d0
      one = 1.d0
      per = period(1) .gt. zero
      if(ix .gt. 0)then
      if(per)then
      call dist2(u,v,x(ix),y(ix),period,d2)
      else
      d2 = (u-x(ix))**2 + (v-y(ix))**2
      endif
      newpt = d2 .gt. eps
      else
      newpt = .true.
      endif
      beta = par(1)
      gamma = par(2)
      r2 = par(3)
      s = par(4)
      if(newpt)then
      c1 = zero
      else
      c1 = dble(aux(ix))
      endif
      c2 = zero
      do23006 j = 1,npts 
      if(j .eq. ix)then
      goto 23006
      endif
      if(ix .gt. 0)then
      if(per)then
      call dist2(x(ix),y(ix),x(j),y(j),period,d2)
      else
      d2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
      endif
      if(d2 .lt. r2)then
      a1 = one
      else
      a1 = zero
      endif
      else
      a1 = zero
      endif
      if(newpt)then
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .lt. r2)then
      a2 = one
      else
      a2 = zero
      endif
      c1 = c1 + a2
      else
      a2 = a1
      endif
      a0 = dble(aux(j))
      c2 = c2 + min(s,a0-a1+a2) - min(s,a0-a1)
23006 continue
23007 continue
      count = min(s,c1) + c2
      if(gamma .lt. eps )then
      if(count .gt. zero)then
      cifval = zero
      else
      cifval = beta
      endif
      else
      cifval = beta*exp(log(gamma)*dble(count))
      endif
      return
      end
