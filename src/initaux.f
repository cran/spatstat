C Output from Public domain Ratfor, version 1.0
      subroutine initaux(nmbr,par,period,x,y,npts,ndisc,aux)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), par(1), period(2)
      integer aux(ndisc,1)
      logical per
      if(nmbr .ne. 8 .and. nmbr .ne. 11)then
      return
      endif
      per = period(1) .gt. 0.d0
      do23002 i=1,npts 
      do23004 k=1,ndisc 
      aux(k,i) = 0
23004 continue
23005 continue
      do23006 j=1,npts 
      if(j.eq.i)then
      goto 23006
      endif
      if(per)then
      call dist2(x(i),y(i),x(j),y(j),period,d2)
      else
      d2 = (x(i)-x(j))**2 + (y(i)-y(j))**2
      endif
      do23012 k = ndisc,1,-1 
      r2 = par(3*k+1)
      if(d2 .lt. r2)then
      aux(k,i) = aux(k,i) + 1
      else
      goto 23013
      endif
23012 continue
23013 continue
23006 continue
23007 continue
23002 continue
23003 continue
      return
      end
