      subroutine straushm(u,v,mrk,ix,x,y,marks,nmarks,n,par,period,
&     cifval)
      implicit double precision(a-h,o-z)
      dimension x(n), y(n), marks(n), par(1), period(2)
      dimension gmma(100), rint(100), rhc(100), k(100)
      logical per
      eps = 2.22d-16
      zero = 0.d0
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
      rint(i) = par(ind)**2
23012 continue
      if(.not.(i.eq.mrk .and. j.gt.mrk))goto 23014
      rint(j) = par(ind)**2
23014 continue
23010 continue
23008 continue
      ind = nmarks + nmarks*(nmarks+1)
      do 23016 i = 1,mrk 
      do 23018 j = i,nmarks 
      ind = ind+1
      if(.not.(j.eq.mrk))goto 23020
      rhc(i) = par(ind)**2
23020 continue
      if(.not.(i.eq.mrk .and. j.gt.mrk))goto 23022
      rhc(j) = par(ind)**2
23022 continue
23018 continue
23016 continue
      do 23024 i = 1,nmarks 
      k(i) = 0
23024 continue
      do 23026 j = 1,n 
      if(.not.(j .eq. ix))goto 23028
      continue
      goto 23029
23028 continue
      if(.not.(per))goto 23030
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23031
23030 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23031 continue
      if(.not.(d2 .lt. rhc(marks(j))))goto 23032
      cifval = zero
      return
23032 continue
      if(.not.(d2 .le. rint(marks(j))))goto 23034
      k(marks(j)) = k(marks(j))+1
23034 continue
23029 continue
23026 continue
      cifval = zero
      do 23036 i=1,nmarks 
      if(.not.(gmma(i) .lt. eps))goto 23038
      if(.not.(k(i).gt.0))goto 23040
      cifval = zero
      return
23040 continue
      goto 23039
23038 continue
      cifval = cifval + log(gmma(i))*dble(k(i))
23039 continue
23036 continue
      cifval = beta*exp(cifval)
      return
      end
