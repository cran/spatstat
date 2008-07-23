C Output from Public domain Ratfor, version 1.0
      subroutine updaux(itype,x,y,u,v,npts,ix,par,period,ndisc,aux)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), par(2), period(2)
      integer aux(ndisc,1)
      logical per
      ndisc = par(2)
      per = period(1) .gt. 0.d0
      nm1 = npts-1
      if(itype .eq. 1)then
      do23002 k=1,ndisc 
      aux(k,npts) = 0
23002 continue
23003 continue
      do23004 j=1,nm1 
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      do23008 k=ndisc,1,-1 
      r2 = par(3*k+1)
      if(d2 .lt. r2)then
      aux(k,j) = aux(k,j)+1
      aux(k,npts) = aux(k,npts)+1
      else
      goto 23009
      endif
23008 continue
23009 continue
23004 continue
23005 continue
      return
      endif
      if(itype .eq. 2)then
      do23014 j = 1,npts 
      if(j.eq.ix)then
      goto 23014
      endif
      if(per)then
      call dist2(x(ix),y(ix),x(j),y(j),period,d2)
      else
      d2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
      endif
      do23020 k=ndisc,1,-1 
      r2 = par(3*k+1)
      if(d2 .ge. r2 .and. j .lt. ix)then
      goto 23020
      endif
      if(d2 .lt. r2)then
      if(j .lt. ix)then
      aux(k,j) = aux(k,j) - 1
      else
      aux(k,j-1) = aux(k,j) - 1
      endif
      else
      if(j.ge.ix)then
      aux(k,j-1) = aux(k,j)
      endif
      endif
23020 continue
23021 continue
23014 continue
23015 continue
      do23030 k = 1,ndisc 
      aux(k,npts) = 0
23030 continue
23031 continue
      return
      endif
      if(itype .eq. 3)then
      do23034 k = 1,ndisc 
      aux(k,ix) = 0
23034 continue
23035 continue
      do23036 j = 1,npts 
      if(j .eq. ix)then
      goto 23036
      endif
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(per)then
      call dist2(x(ix),y(ix),x(j),y(j),period,d2i)
      else
      d2i = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
      endif
      do23044 k = 1,ndisc 
      r2 = par(3*k+1)
      if(d2 .ge. r2 .and. d2i .ge. r2)then
      goto 23044
      endif
      if(d2 .lt. r2)then
      aux(k,ix) = aux(k,ix) + 1
      if(d2i .ge. r2)then
      aux(k,j) = aux(k,j) + 1
      endif
      else
      if(d2i .lt. r2)then
      aux(k,j) = aux(k,j) - 1
      endif
      endif
23044 continue
23045 continue
23036 continue
23037 continue
      return
      endif
      call fexit("Argument itype to updaux must be 1, 2, or 3; bailing o
     *ut.\n")
      end
