C Output from Public domain Ratfor, version 1.0
      subroutine initaux(nmbr,par,period,x,y,npts,aux)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), par(3), period(2)
      integer aux(1)
      logical per
      if(nmbr .ne. 8)then
      return
      endif
      per = period(1) .gt. 0.d0
      r2 = par(3)**2
      do23002 i = 1, npts 
      aux(i) = 0
      do23004 j = 1,npts 
      if(j.eq.i)then
      goto 23004
      endif
      if(per)then
      call dist2(x(i),y(i),x(j),y(j),period,d2)
      else
      d2 = (x(i)-x(j))**2 + (y(i)-y(j))**2
      endif
      if(d2 .lt. r2)then
      aux(i) = aux(i) + 1
      endif
23004 continue
23005 continue
23002 continue
23003 continue
      return
      end
