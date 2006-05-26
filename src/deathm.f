C Output from Public domain Ratfor, version 1.0
      subroutine deathm(x,y,marks,npts,ix)
      implicit double precision(a-h,o-z)
      dimension x(1), y(1)
      dimension marks(1)
      npts = npts-1
      if(ix .gt. npts)then
      return
      endif
      do23002 j = ix,npts 
      x(j) = x(j+1)
      y(j) = y(j+1)
      marks(j) = marks(j+1)
23002 continue
23003 continue
      return
      end
