      subroutine arand(ix,iy,iz,rand)
      implicit double precision(a-h,o-z)
      ix = 171 * (ix - 177*(ix/177)) - 2 * (ix/177)
      iy = 172 * (iy - 176*(iy/176)) - 35 * (iy/176)
      iz = 170 * (iz - 178*(iz/178)) - 63 * (iz/178)
      if(.not.(ix .lt. 0))goto 23000
      ix = ix + 30269
23000 continue
      if(.not.(iy .lt. 0))goto 23002
      iy = iy + 30307
23002 continue
      if(.not.(iz .lt. 0))goto 23004
      iz = iz + 30323
23004 continue
      term = ix/30269.d0 + iy/30307.d0 + iz/30323.d0
      iterm = term
      if(.not.(iterm .gt. term))goto 23006
      iterm = iterm - 1
23006 continue
      rand = term - iterm
      if(.not.(rand .lt. 1.d0))goto 23008
      continue
      goto 23009
23008 continue
      rand = 0.999999d0
23009 continue
      return
      end
