C Output from Public domain Ratfor, version 1.0
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
      do23000 j = 1,n 
      if(j .eq. ix)then
      continue
      else
      if(per)then
      call dist2(u,v,x(j),y(j),period,d2)
      else
      d2 = (u-x(j))**2 + (v-y(j))**2
      endif
      sx = sx + exp(oneomk*log(d2))
      endif
23000 continue
23001 continue
      cifval = beta*exp(-(sigma**twook)*sx)
      return
      end
