      subroutine aru(n,a,b,iseed,rrr)
      implicit double precision(a-h,o-z)
      dimension iseed(3), rrr(n)
      w = b-a
      do 23000 i = 1,n 
      call arand(iseed(1),iseed(2),iseed(3),rv)
      rrr(i) = a + w*rv
23000 continue
      return
      end
