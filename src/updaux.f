C Output from Public domain Ratfor, version 1.0
      subroutine updaux(itype,x,y,u,v,npts,ix,par,period,aux)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), par(4), period(2)
      integer aux(1)
      logical per
      r2 = par(3)**2
      per = period(1) .gt. 0.d0
      if(itype .eq. 1)then
      nm1 = npts-1
      aux(npts) = 0
      do23002 j = 1,nm1 
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      if(d2 .lt. r2)then
      aux(j) = aux(j)+1
      aux(npts) = aux(npts)+1
      endif
23002 continue
23003 continue
      return
      endif
      if(itype .eq. 2)then
      do23010 j = 1,npts 
      if(j.eq.ix)then
      goto 23010
      endif
      if(per)then
      call dist2(x(ix),y(ix),x(j),y(j),period,d2)
      else
      d2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
      endif
      if(d2 .lt. r2)then
      if(j .lt. ix)then
      aux(j) = aux(j) - 1
      else
      aux(j-1) = aux(j) - 1
      endif
      else
      if(j.ge.ix)then
      aux(j-1) = aux(j)
      endif
      endif
23010 continue
23011 continue
      aux(npts) = 0
      return
      endif
      if(itype .eq. 3)then
      aux(ix) = 0
      do23024 j = 1, npts 
      if(j .eq. ix)then
      goto 23024
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
      if(d2 .lt. r2)then
      aux(ix) = aux(ix) + 1
      if(d2i .ge. r2)then
      aux(j) = aux(j) + 1
      endif
      else
      if(d2i .lt. r2)then
      aux(j) = aux(j) - 1
      endif
      endif
23024 continue
23025 continue
      return
      endif
      call fexit("Argument itype to updaux must be 1, 2, or 3; bailing o
     *ut.")
      end
