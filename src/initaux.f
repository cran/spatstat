      subroutine initaux(nmbr,par,period,x,y,npts,aux)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), par(4), period(2)
      integer aux(1)
      logical per
      if(.not.(nmbr .ne. 8))goto 23000
      return
23000 continue
      per = period(1) .gt. 0.d0
      r2 = par(3)**2
      do 23002 i = 1, npts 
      aux(i) = 0
      do 23004 j = 1,npts 
      if(.not.(j.eq.i))goto 23006
      goto 23004
23006 continue
      if(.not.(per))goto 23008
      call dist2(x(i),y(i),x(j),y(j),period,d2)
      goto 23009
23008 continue
      d2 = (x(i)-x(j))**2 + (y(i)-y(j))**2
23009 continue
      if(.not.(d2 .lt. r2))goto 23010
      aux(i) = aux(i) + 1
23010 continue
23004 continue
23002 continue
      return
      end
