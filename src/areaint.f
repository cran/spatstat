C Output from Public domain Ratfor, version 1.0
      subroutine areaint(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      dimension par(3), x(1), y(1), period(2)
      logical per
      logical covered
      logical neigh(n)
      eps = 2.22d-16
      zero = 0.d0
      ngrid = 64
      per = period(1) .gt. zero
      beta = par(1)
      eta = par(2)
      r = par(3)
      if(n.eq.0)then
      cifval = beta
      return
      endif
      r2 = r * r
      dx = (2 * r)/dble(ngrid)
      dy = dx
      range2 = 4 * r2
      if(.not. per)then
      ixm1 = ix - 1
      ixp1 = max(1, ix + 1)
      if(ixm1 .gt. 0)then
      do23006 j = 1,ixm1 
      a = range2 - (u - x(j))**2
      neigh(j) = .false.
      if(a .gt. 0)then
      a = a - (v - y(j))**2
      if(a .gt. 0)then
      neigh(j) = .true.
      endif
      endif
23006 continue
23007 continue
      endif
      if(ixp1 .le. n)then
      do23014 j = ixp1,n 
      a = range2 - (u - x(j))**2
      neigh(j) = .false.
      if(a .gt. 0)then
      a = a - (v - y(j))**2
      if(a .gt. 0)then
      neigh(j) = .true.
      endif
      endif
23014 continue
23015 continue
      endif
      kount = 0
      kdisc = 0
      xgrid0 = u - r - dx/2
      do23020 kx=1,ngrid 
      xgrid = xgrid0 + dble(kx) * dx
      my = int(dsqrt(r2 - (u - xgrid)**2)/dy)
      if(my .ge. 0)then
      do23024 ky=-my,my 
      ygrid = v + dble(ky) * dy
      kdisc = kdisc + 1
      covered = .false.
      if(ixm1 .gt. 0)then
      do23028 j = 1,ixm1 
      if(neigh(j))then
      a = r2 - (xgrid - x(j))**2
      if(a .gt. 0)then
      a = a - (ygrid - y(j))**2
      if(a .gt. 0)then
      covered = .true.
      goto 42
      endif
      endif
      endif
23028 continue
23029 continue
      endif
      if(ixp1 .le. n)then
      do23038 j = ixp1,n 
      if(neigh(j))then
      a = r2 - (xgrid - x(j))**2
      if(a .gt. 0)then
      a = a - (ygrid - y(j))**2
      if(a .gt. 0)then
      covered = .true.
      goto 42
      endif
      endif
      endif
23038 continue
23039 continue
      endif
      if(.not. covered)then
      kount = kount + 1
      endif
42    continue
23024 continue
23025 continue
      endif
23020 continue
23021 continue
      else
      ixm1 = ix - 1
      ixp1 = max(1, ix + 1)
      if(ixm1 .gt. 0)then
      do23050 j = 1,ixm1 
      call dist2(u,v,x(j),y(j),period,d2)
      neigh(j) = (d2 .lt. range2)
23050 continue
23051 continue
      endif
      if(ixp1 .le. n)then
      do23054 j = ixp1,n 
      call dist2(u,v,x(j),y(j),period,d2)
      neigh(j) = (d2 .lt. range2)
23054 continue
23055 continue
      endif
      kount = 0
      kdisc = 0
      xgrid0 = u - r - dx
      ygrid0 = v - r - dy
      do23056 kx=1,ngrid 
      xgrid = xgrid0 + dble(kx) * dx
      do23058 ky=1,ngrid 
      ygrid = ygrid0 + dble(ky) * dy
      call dist2(u,v,xgrid,ygrid,period,d2)
      if(d2 .lt. r2)then
      kdisc = kdisc + 1
      covered = .false.
      if(ixm1 .gt. 0)then
      do23064 j = 1,ixm1 
      if(neigh(j))then
      call dist2(xgrid,ygrid,x(j),y(j),period,d2)
      if(d2 .lt. r2)then
      covered = .true.
      goto 4242
      endif
      endif
23064 continue
23065 continue
      endif
      if(ixp1 .le. n)then
      do23072 j = ixp1,n 
      if(neigh(j))then
      call dist2(xgrid,ygrid,x(j),y(j),period,d2)
      if(d2 .lt. r2)then
      covered = .true.
      goto 4242
      endif
      endif
23072 continue
23073 continue
      endif
      if(.not. covered)then
      kount = kount + 1
      endif
4242  continue
      endif
23058 continue
23059 continue
23056 continue
23057 continue
      endif
      if(eta .gt. eps)then
      covfrac = dble(kdisc - kount) / dble(kdisc)
      cifval = beta * exp(log(eta) * covfrac)
      else
      if(kount .eq. kdisc)then
      cifval = beta
      else
      cifval = zero
      endif
      endif
      return
      end
