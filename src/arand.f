C Output from Public domain Ratfor, version 1.0
      subroutine arand(ix,iy,iz,rand)
      implicit double precision(a-h,o-z)
      ix = 171 * (ix - 177*(ix/177)) - 2 * (ix/177)
      iy = 172 * (iy - 176*(iy/176)) - 35 * (iy/176)
      iz = 170 * (iz - 178*(iz/178)) - 63 * (iz/178)
      if(ix .lt. 0)then
      ix = ix + 30269
      endif
      if(iy .lt. 0)then
      iy = iy + 30307
      endif
      if(iz .lt. 0)then
      iz = iz + 30323
      endif
      term = ix/30269.d0 + iy/30307.d0 + iz/30323.d0
      iterm = term
      if(iterm .gt. term)then
      iterm = iterm - 1
      endif
      rand = term - iterm
      if(rand .lt. 1.d0)then
      continue
      else
      rand = 0.999999d0
      endif
      return
      end
