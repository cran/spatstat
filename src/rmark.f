      subroutine rmark(n,ntypes,ptypes,iseed,marks)
      implicit double precision(a-h,o-z)
      dimension iseed(3), ptypes(1), marks(1)
      do 23000 i = 1,n 
      call arand(iseed(1),iseed(2),iseed(3),rv)
      psum = 0.d0
      do 23002 j = 1,ntypes 
      psum = psum+ptypes(j)
      if(.not.(rv .lt. psum))goto 23004
      marks(i) = j
      goto 23003
23004 continue
23002 continue
23003 continue
23000 continue
      return
      end
