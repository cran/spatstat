      subroutine death(x,y,marks,marked,npts,ix)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), marks(1)
      logical marked
      npts = npts-1
      if(.not.(ix .gt. npts))goto 23000
      return
23000 continue
      do 23002 j = ix,npts 
      x(j) = x(j+1)
      y(j) = y(j+1)
      if(.not.(marked))goto 23004
      marks(j) = marks(j+1)
23004 continue
23002 continue
      return
      end
