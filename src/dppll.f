      subroutine dppll(x,y,l1,l2,l3,l4,np,nl,eps,mint,rslt,xmin,jmin)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), rslt(np,1), xmin(1), jmin(1)
      double precision l1(1), l2(1), l3(1), l4(1)
      one = 1.d0
      zero = 0.d0
      do 23000 j = 1,nl 
      dx = l3(j) - l1(j)
      dy = l4(j) - l2(j)
      alen = sqrt(dx**2 + dy**2)
      if(.not.(alen .gt. eps))goto 23002
      co = dx/alen
      si = dy/alen
23002 continue
      do 23004 i = 1, np 
      xpx1 = x(i) - l1(j)
      ypy1 = y(i) - l2(j)
      xpx2 = x(i) - l3(j)
      ypy2 = y(i) - l4(j)
      d1 = xpx1**2 + ypy1**2
      d2 = xpx2**2 + ypy2**2
      dd = min(d1,d2)
      if(.not.(alen .gt. eps))goto 23006
      xpr = xpx1*co + ypy1*si
      if(.not.(xpr .lt. zero .or. xpr .gt. alen))goto 23008
      d3 = -one
      goto 23009
23008 continue
      ypr = - xpx1*si + ypy1*co
      d3 = ypr**2
23009 continue
      goto 23007
23006 continue
      d3 = -one
23007 continue
      if(.not.(d3 .gt. zero))goto 23010
      dd = min(dd,d3)
23010 continue
      sd =sqrt(dd)
      rslt(i,j) = sd
      if(.not.(mint.gt.0))goto 23012
      if(.not.(sd .lt. xmin(i)))goto 23014
      xmin(i) = sd
      if(.not.(mint.gt.1))goto 23016
      jmin(i) = j
23016 continue
23014 continue
23012 continue
23004 continue
23000 continue
      return
      end
