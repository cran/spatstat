C Output from Public domain Ratfor, version 1.0
      subroutine straussm(u,v,mrk,ix,x,y,marks,n,ntypes,par,period,cifva
     *l)
      implicit double precision(a-h,o-z)
      dimension x(n), y(n), marks(n), par(1), period(2)
      dimension gmma(100), rad(100), k(100)
      logical per
      eps = 2.22d-16
      zero = 0.d0
      per = period(1) .gt. zero
      beta = par(mrk)
      ind = ntypes
      do23000 i = 1,mrk 
      do23002 j = i,ntypes 
      ind = ind+1
      if(i.le.mrk .and. j.eq.mrk)then
      gmma(i) = par(ind)
      endif
      if(i.eq.mrk .and. j.gt.mrk)then
      gmma(j) = par(ind)
      endif
23002 continue
23003 continue
23000 continue
23001 continue
      ind = ntypes + ntypes*(ntypes+1)/2
      do23008 i = 1,mrk 
      do23010 j = i,ntypes 
      ind = ind+1
      if(i.le.mrk .and. j.eq.mrk)then
      rad(i) = par(ind)
      endif
      if(i.eq.mrk .and. j.gt.mrk)then
      rad(j) = par(ind)
      endif
23010 continue
23011 continue
23008 continue
23009 continue
      do23016 i = 1,ntypes 
      k(i) = 0
23016 continue
23017 continue
      do23018 j = 1,n 
      if(j .eq. ix)then
      continue
      else
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .le. rad(marks(j)))then
      k(marks(j)) = k(marks(j))+1
      endif
      endif
23018 continue
23019 continue
      cifval = zero
      do23026 i=1,ntypes 
      if(gmma(i) .lt. eps)then
      if(k(i).gt.0)then
      cifval = zero
      return
      endif
      else
      cifval = cifval + log(gmma(i))*dble(k(i))
      endif
23026 continue
23027 continue
      cifval = beta*exp(cifval)
      return
      end
