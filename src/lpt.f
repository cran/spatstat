C Output from Public domain Ratfor, version 1.0
      subroutine lpt(u,v,tpar,tval)
      implicit double precision(a-h,o-z)
      dimension tpar(1)
      poly = 0.d0
      nd = nint(tpar(1))
      i = 1
      do23000 j = 1,nd 
      do23002 k = 0,j 
      i = i + 1
      poly = poly + tpar(i)*u**(j-k)*v**k
23002 continue
23003 continue
23000 continue
23001 continue
      tval = exp(poly)
      return
      end
