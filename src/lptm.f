C Output from Public domain Ratfor, version 1.0
      subroutine lptm(u,v,mark,tpar,tval)
      implicit double precision(a-h,o-z)
      dimension tpar(1)
      nmarks = nint(tpar(1))
      nstart = 1 + nmarks
      if(mark .gt. 1)then
      do23002 i = 2,mark 
      nd = nint(tpar(i))
      nc = nd + nd*(nd+1)/2
      nstart = nstart + nc
23002 continue
23003 continue
      endif
      poly = 0.d0
      i = nstart
      nd = nint(tpar(1+mark))
      if(nd .gt. 0)then
      do23006 j = 1,nd 
      do23008 k = 0,j 
      i = i + 1
      poly = poly + tpar(i)*u**(j-k)*v**k
23008 continue
23009 continue
23006 continue
23007 continue
      endif
      tval = exp(poly)
      return
      end
