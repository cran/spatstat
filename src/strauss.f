      subroutine strauss(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      dimension par(3), x(1), y(1), period(2)
      logical per
      eps = 2.22d-16
      per = period(1) .gt. 0.d0
      beta = par(1)
      gamma = par(2)
      r = par(3)**2
      kount = 0
      do 23000 j = 1,n 
      if(.not.(j .eq. ix))goto 23002
      continue
      goto 23003
23002 continue
      if(.not.(per))goto 23004
      call dist2(u,v,x(j),y(j),period,d2)
      goto 23005
23004 continue
      d2 = (u-x(j))**2 + (v-y(j))**2
23005 continue
      if(.not.(d2 .lt. r))goto 23006
      kount = kount+1
23006 continue
23003 continue
23000 continue
      if(.not.(gamma .lt. eps ))goto 23008
      if(.not.(kount .gt. 0))goto 23010
      cifval = 0.d0
      goto 23011
23010 continue
      cifval = beta
23011 continue
      goto 23009
23008 continue
      cifval = beta*exp(log(gamma)*dble(kount))
23009 continue
      return
      end
