C Output from Public domain Ratfor, version 1.0
      subroutine dist2(u,v,x,y,period,d2)
      implicit double precision(a-h,o-z)
      dimension period(2)
      wide = period(1)
      high = period(2)
      a1 = (u-x)**2
      a2 = (u-x+wide)**2
      a3 = (u-x-wide)**2
      b1 = (v-y)**2
      b2 = (v-y+high)**2
      b3 = (v-y-high)**2
      amin = min(a1,a2,a3)
      bmin = min(b1,b2,b3)
      d2 = amin + bmin
      return
      end
