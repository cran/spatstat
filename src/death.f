C Output from Public domain Ratfor, version 1.0
      subroutine death(x,y,marks,marked,npts,ix)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1), marks(1)
      logical marked
      npts = npts-1
      if(ix .gt. npts)then
      return
      endif
      do23002 j = ix,npts 
      x(j) = x(j+1)
      y(j) = y(j+1)
      if(marked)then
      marks(j) = marks(j+1)
      endif
23002 continue
23003 continue
      return
      end
