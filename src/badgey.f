C Output from Public domain Ratfor, version 1.0
      subroutine badgey(u,v,ix,x,y,npts,par,period,cifval,ndisc,aux)
      implicit double precision (a-h,o-z)
      integer aux(ndisc,1)
      dimension x(1), y(1), marks(1), par(1), period(2)
      dimension w(ndisc), tee(ndisc)
      logical newpt, per
      beta = par(1)
      if(npts.eq.0)then
      cifval = beta
      return
      endif
      eps = 2.22e-16
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
      do23006 i=1,ndisc 
      w(i) = zero
23006 continue
23007 continue
      if(newpt)then
      do23010 i=1,ndisc 
      tee(i) = zero
23010 continue
23011 continue
      do23012 j=1,npts 
      if(j .eq. ix)then
      goto 23012
      endif
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      do23018 i=ndisc,1,-1 
      r2 = par(3*i+1)
      s = par(3*i+2)
      if(d2 .lt. r2)then
      tee(i) = tee(i) + one
      a = aux(i,j)
      if(ix.gt.0)then
      if(per)then
      call dist2(x(ix),y(ix), x(j),y(j),period,dd2)
      else
      dd2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
      endif
      if(dd2 .lt. r2)then
      a = a - one
      endif
      endif
      b = a + one
      if(a .lt. s .and. s .lt. b)then
      w(i) = w(i) + s - a
      else
      if(s .ge. b)then
      w(i) = w(i) + one
      endif
      endif
      else
      goto 23019
      endif
23018 continue
23019 continue
23012 continue
23013 continue
      else
      do23032 i=1,ndisc 
      tee(i) = aux(i,ix)
23032 continue
23033 continue
      do23034 j=1,npts 
      if(j .eq. ix)then
      goto 23034
      endif
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      do23040 i=ndisc,1,-1 
      r2 = par(3*i+1)
      s = par(3*i+2)
      if(d2 .lt. r2)then
      a = aux(i,j) - 1
      b = a + one
      if(a .lt. s .and. s .lt. b)then
      w(i) = w(i) + s - a
      else
      if(s .ge. b)then
      w(i) = w(i) + one
      endif
      endif
      else
      goto 23041
      endif
23040 continue
23041 continue
23034 continue
23035 continue
      endif
      some = zero
      do23048 i=1,ndisc 
      s = par(3*i+2)
      w(i) = w(i) + min(tee(i),s)
      gamma = par(3*i)
      if(gamma .lt. eps )then
      if(w(i) .gt. zero)then
      cifval = zero
      return
      endif
      else
      some = some + log(gamma)*w(i)
      endif
23048 continue
23049 continue
      cifval = beta*exp(some)
      return
      end
