      subroutine sftcr(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      double precision kappa
      dimension par(3), x(n), y(n), period(2)
      logical per
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      beta = par(1)
      sigma = par(2)
      kappa = par(3)
      oneomk = -one/kappa
      twook = two/kappa
      per = period(1) .gt. 0.d0
      sx = zero
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
      sx = sx + exp(oneomk*log(d2))
23003 continue
23000 continue
      cifval = beta*exp(-(sigma**twook)*sx)
      return
      end
