C Output from Public domain Ratfor, version 1.0
      subroutine aru(n,a,b,iseed,rrr)
      implicit double precision(a-h,o-z)
      dimension iseed(3), rrr(n)
      w = b-a
      do23000 i = 1,n 
      call arand(iseed(1),iseed(2),iseed(3),rv)
      rrr(i) = a + w*rv
23000 continue
23001 continue
      return
      end
