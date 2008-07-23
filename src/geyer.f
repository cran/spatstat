C Output from Public domain Ratfor, version 1.0
      subroutine geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
      implicit double precision(a-h,o-z)
      dimension par(5), x(1), y(1), period(2)
      integer aux(1)
      logical per, newpt
      eps = 2.22d-16
      zero = 0.d0
      one = 1.d0
      per = period(1) .gt. zero
      beta = par(1)
      gamma = par(3)
      r2 = par(4)
      s = par(5)
      if(npts.eq.0)then
      cifval = beta
      return
      endif
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
      w = zero
      if(newpt)then
      tee = zero
      do23008 j = 1,npts 
      if(j .eq. ix)then
      goto 23008
      endif
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .lt. r2)then
      tee = tee + one
      a = aux(j)
      if(ix .gt. 0)then
      if(per)then
      call dist2(x(ix),y(ix),x(j),y(j),period,dd2)
      else
      dd2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
      endif
      if(dd2 .lt. r2)then
      a = a - one
      endif
      endif
      b = a + one
      if(a .lt. s .and. s .lt. b)then
      w = w + s - a
      else
      if(s .ge. b)then
      w = w + one
      endif
      endif
      endif
23008 continue
23009 continue
      else
      tee = aux(ix)
      do23026 j = 1,npts 
      if(j .eq. ix)then
      goto 23026
      endif
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .lt. r2)then
      a = aux(j) - 1
      b = a + one
      if(a .lt. s .and. s .lt. b)then
      w = w + s - a
      else
      if(s .ge. b)then
      w = w + one
      endif
      endif
      endif
23026 continue
23027 continue
      endif
      w = w + min(tee,s)
      count = min(s,c1) + c2
      if(gamma .lt. eps )then
      if(w .gt. zero)then
      cifval = zero
      else
      cifval = beta
      endif
      else
      cifval = beta*exp(log(gamma)*w)
      endif
      return
      end
