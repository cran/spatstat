      subroutine diggra(u,v,ix,x,y,n,par,period,cifval)
      implicit double precision(a-h,o-z)
      double precision kappa
      dimension par(4), x(n), y(n), period(2)
      logical per
      zero = 0.d0
      per = period(1) .gt. 0.d0
      beta = par(1)
      kappa = par(2)
      delta = par(3)
      rho = par(4)
      d2 = delta**2
      r2 = rho**2
      a = log(rho - delta)
      soglum = zero
      do 23000 j = 1,n 
      if(.not.(j .eq. ix))goto 23002
      sincr = zero
      goto 23003
23002 continue
      if(.not.(per))goto 23004
      call dist2(u,v,x(j),y(j),period,t2)
      goto 23005
23004 continue
      t2 = (u-x(j))**2 + (v-y(j))**2
23005 continue
      if(.not.(t2 .lt. d2))goto 23006
      cifval = zero
      return
23006 continue
      if(.not.(t2 .lt. r2))goto 23008
      sincr = log(sqrt(t2) - delta) - a
      goto 23009
23008 continue
      sincr = zero
23009 continue
23003 continue
      soglum = soglum + sincr
23000 continue
      cifval = beta*exp(kappa*soglum)
      return
      end
