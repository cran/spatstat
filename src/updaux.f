      subroutine updaux(itype,x,y,u,v,npts,ix,par,period,aux)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), par(4), period(2)
      integer aux(1)
      logical per
      r2 = par(3)**2
      per = period(1) .gt. 0.d0
      if(.not.(itype .eq. 1))goto 23000
      nm1 = npts-1
      aux(npts) = 0
      do 23002 j = 1,nm1 
      if(.not.(per))goto 23004
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23005
23004 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23005 continue
      if(.not.(d2 .lt. r2))goto 23006
      aux(j) = aux(j)+1
      aux(npts) = aux(npts)+1
23006 continue
23002 continue
      return
23000 continue
      if(.not.(itype .eq. 2))goto 23008
      do 23010 j = 1,npts 
      if(.not.(j.eq.ix))goto 23012
      goto 23010
23012 continue
      if(.not.(per))goto 23014
      call dist2(x(ix),y(ix),x(j),y(j),period,d2)
      goto 23015
23014 continue
      d2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
23015 continue
      if(.not.(d2 .lt. r2))goto 23016
      if(.not.(j .lt. ix))goto 23018
      aux(j) = aux(j) - 1
      goto 23019
23018 continue
      aux(j-1) = aux(j) - 1
23019 continue
      goto 23017
23016 continue
      if(.not.(j.ge.ix))goto 23020
      aux(j-1) = aux(j)
23020 continue
23017 continue
23010 continue
      aux(npts) = 0
      return
23008 continue
      if(.not.(itype .eq. 3))goto 23022
      aux(ix) = 0
      do 23024 j = 1, npts 
      if(.not.(j .eq. ix))goto 23026
      goto 23024
23026 continue
      if(.not.(per))goto 23028
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23029
23028 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23029 continue
      if(.not.(per))goto 23030
      call dist2(x(ix),y(ix),x(j),y(j),period,d2i)
      goto 23031
23030 continue
      d2i = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
23031 continue
      if(.not.(d2 .lt. r2))goto 23032
      aux(ix) = aux(ix) + 1
      if(.not.(d2i .ge. r2))goto 23034
      aux(j) = aux(j) + 1
23034 continue
      goto 23033
23032 continue
      if(.not.(d2i .lt. r2))goto 23036
      aux(j) = aux(j) - 1
23036 continue
23033 continue
23024 continue
      return
23022 continue
      call fexit(
&     "Argument itype to updaux must be 1, 2, or 3; bailing out.")
      end
