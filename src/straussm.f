      subroutine straussm(u,v,mrk,ix,x,y,marks,n,nmarks,par,period,
&     cifval)
      implicit double precision(a-h,o-z)
      dimension x(n), y(n), marks(n), par(1), period(2)
      dimension gmma(100), rad(100), k(100)
      logical per
      eps = 2.22d-16
      per = period(1) .gt. 0.d0
      beta = par(mrk)
      ind = nmarks
      do 23000 i = 1,mrk 
      do 23002 j = i,nmarks 
      ind = ind+1
      if(.not.(i.le.mrk .and. j.eq.mrk))goto 23004
      gmma(i) = par(ind)
23004 continue
      if(.not.(i.eq.mrk .and. j.gt.mrk))goto 23006
      gmma(j) = par(ind)
23006 continue
23002 continue
23000 continue
      ind = nmarks + nmarks*(nmarks+1)/2
      do 23008 i = 1,mrk 
      do 23010 j = i,nmarks 
      ind = ind+1
      if(.not.(i.le.mrk .and. j.eq.mrk))goto 23012
      rad(i) = par(ind)**2
23012 continue
      if(.not.(i.eq.mrk .and. j.gt.mrk))goto 23014
      rad(j) = par(ind)**2
23014 continue
23010 continue
23008 continue
      do 23016 i = 1,nmarks 
      k(i) = 0
23016 continue
      do 23018 j = 1,n 
      if(.not.(j .eq. ix))goto 23020
      continue
      goto 23021
23020 continue
      if(.not.(per))goto 23022
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23023
23022 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23023 continue
      if(.not.(d2 .le. rad(marks(j))))goto 23024
      k(marks(j)) = k(marks(j))+1
23024 continue
23021 continue
23018 continue
      cifval = zero
      do 23026 i=1,nmarks 
      if(.not.(gmma(i) .lt. eps))goto 23028
      if(.not.(k(i).gt.0))goto 23030
      cifval = 0.d0
      return
23030 continue
      goto 23029
23028 continue
      cifval = cifval + log(gmma(i))*dble(k(i))
23029 continue
23026 continue
      cifval = beta*exp(cifval)
      return
      end
