      subroutine geyer(u,v,ix,x,y,npts,par,period,cifval,aux)
      implicit double precision(a-h,o-z)
      dimension par(4), x(1), y(1), period(2)
      integer aux(1)
      logical per, newpt
      eps = 2.22d-16
      zero = 0.d0
      one = 1.d0
      per = period(1) .gt. 0.d0
      if(.not.(ix .gt. 0))goto 23000
      if(.not.(per))goto 23002
      call dist2(u,v,x(ix),y(ix),period,d2)
      goto 23003
23002 continue
      d2 = (u-x(ix))**2 + (v-y(ix))**2
23003 continue
      newpt = d2 .gt. eps
      goto 23001
23000 continue
      newpt = .true.
23001 continue
      beta = par(1)
      gamma = par(2)
      r2 = par(3)**2
      s = par(4)
      if(.not.(newpt))goto 23004
      c1 = zero
      goto 23005
23004 continue
      c1 = dble(aux(ix))
23005 continue
      c2 = zero
      do 23006 j = 1,npts 
      if(.not.(j .eq. ix))goto 23008
      goto 23006
23008 continue
      if(.not.(ix .gt. 0))goto 23010
      if(.not.(per))goto 23012
      call dist2(x(ix),y(ix),x(j),y(j),period,d2)
      goto 23013
23012 continue
      d2 = (x(ix)-x(j))**2 + (y(ix)-y(j))**2
23013 continue
      if(.not.(d2 .lt. r2))goto 23014
      a1 = one
      goto 23015
23014 continue
      a1 = zero
23015 continue
      goto 23011
23010 continue
      a1 = zero
23011 continue
      if(.not.(newpt))goto 23016
      if(.not.(per))goto 23018
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23019
23018 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23019 continue
      if(.not.(d2 .lt. r2))goto 23020
      a2 = one
      goto 23021
23020 continue
      a2 = zero
23021 continue
      c1 = c1 + a2
      goto 23017
23016 continue
      a2 = a1
23017 continue
      a0 = dble(aux(j))
      c2 = c2 + min(s,a0-a1+a2) - min(s,a0-a1)
23006 continue
      count = min(s,c1) + c2
      if(.not.(gamma .lt. eps ))goto 23022
      if(.not.(count .gt. zero))goto 23024
      cifval = zero
      goto 23025
23024 continue
      cifval = beta
23025 continue
      goto 23023
23022 continue
      cifval = beta*exp(log(gamma)*dble(count))
23023 continue
      return
      end
