      subroutine lptm(u,v,mark,tpar,tval)
      implicit double precision(a-h,o-z)
      dimension tpar(1)
      nmarks = nint(tpar(1))
      nstart = 1 + nmarks
      if(.not.(mark .gt. 1))goto 23000
      do 23002 i = 2,mark 
      nd = nint(tpar(i))
      nc = nd + nd*(nd+1)/2
      nstart = nstart + nc
23002 continue
23000 continue
      poly = 0.d0
      i = nstart
      nd = nint(tpar(1+mark))
      if(.not.(nd .gt. 0))goto 23004
      do 23006 j = 1,nd 
      do 23008 k = 0,j 
      i = i + 1
      poly = poly + tpar(i)*u**(j-k)*v**k
23008 continue
23006 continue
23004 continue
      tval = exp(poly)
      return
      end
