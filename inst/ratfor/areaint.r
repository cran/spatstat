# $Id: areaint.r,v 1.8 2008/02/11 16:09:18 adrian Exp adrian $
subroutine areaint(u,v,ix,x,y,n,par,period,cifval)
#
# Conditional intensity function for area-interaction process
#
# cif = beta * eta^(-area/pi r^2) = beta * etainv^(area/pi r^2)
#
implicit double precision(a-h,o-z)
dimension par(3), x(1), y(1), period(2)
logical per
logical covered

# dynamically allocated
logical neigh(n)

##
eps  = 2.22d-16 # Essentially .Machine$double.eps 
zero = 0.d0
pi   = 3.14159265359d0
ngrid = 64 

per  = period(1) > zero

beta   = par(1)
etainv = par(2)
r      = par(3)

if(n==0) {
  # uncovered area is pi r^2
  cifval = beta * etainv
  return
}

r2     = r * r
dx     = (2 * r)/dble(ngrid)  # sic
dy     = dx
range2 = 4 * r2     # square of interaction distance 

if(.not. per) {
  # Euclidean distance
  # First identify which data points are neighbours of (u,v)
  ixm1 = ix - 1
  ixp1 = max(1, ix + 1)
  if(ixm1 > 0) {
    do j = 1,ixm1 {
      a = range2 - (u - x(j))**2
      neigh(j) = .false.
      if(a > 0) {
        a = a - (v - y(j))**2
        if(a > 0)
          neigh(j) = .true.
      }
    }
  }
  if(ixp1 <= n) {
    do j = ixp1,n {
      a = range2 - (u - x(j))**2
      neigh(j) = .false.
      if(a > 0) {
        a = a - (v - y(j))**2
        if(a > 0)
          neigh(j) = .true.
      }
    }
  }
  # scan a grid of points centred at (u,v)
  kount = 0
  xgrid0 = u - r - dx/2
  do kx=1,ngrid {
    xgrid = xgrid0 + dble(kx) * dx
    my = int(dsqrt(r2 - (u - xgrid)**2)/dx)
    if(my >= 0) {
      do ky=-my,my {
        ygrid = v + dble(ky) * dy
        # Grid point is inside disc of radius r centred at (u,v)
        # Loop through all data points to determine
        # whether the grid point is covered by another disc
        covered = .false.
        if(ixm1 > 0) {
          do j = 1,ixm1 {
            if(neigh(j)) {
              a = r2 - (xgrid - x(j))**2
              if(a > 0) {
                a = a - (ygrid - y(j))**2
                if(a > 0) {
                  # point j covers grid point
                  covered = .true.
                  goto 42
                }
              }
            }
          }
        }
        if(ixp1 <= n) {
          do j = ixp1,n {
            if(neigh(j)) {
              a = r2 - (xgrid - x(j))**2
              if(a > 0) {
                a = a - (ygrid - y(j))**2
                if(a > 0) {
                  # point j covers grid point
                  covered = .true.
                  goto 42
                }
              }
            }
          }
        }
      }
    }
  }
  if(.not. covered)
    kount = kount + 1
42 continue
} else {
  # periodic distance
  # First identify which data points are neighbours of (u,v)
  ixm1 = ix - 1
  ixp1 = max(1, ix + 1)
  if(ixm1 > 0) {
    do j = 1,ixm1 {
      call dist2(u,v,x(j),y(j),period,d2)
      neigh(j) = (d2 < range2)
    }
  }
  if(ixp1 <= n) {
    do j = ixp1,n {
      call dist2(u,v,x(j),y(j),period,d2)
      neigh(j) = (d2 < range2)
    }
  }
  # scan a grid of ngrid * ngrid points centred at (u,v)
  kount = 0
  xgrid0 = u - r - dx
  ygrid0 = v - r - dy
  do kx=1,ngrid {
    xgrid = xgrid0 + dble(kx) * dx
    do ky=1,ngrid {
      ygrid = ygrid0 + dble(ky) * dy
      # Determine whether grid point is inside disc
      call dist2(u,v,xgrid,ygrid,period,d2)
      if(d2 < r2) {
        # Grid point is inside disc of radius r centred at (u,v)
        # Loop through all data points to determine
        # whether the grid point is covered by another disc
        covered = .false.
        if(ixm1 > 0) {
          do j = 1,ixm1 {
            if(neigh(j)) {
              call dist2(xgrid,ygrid,x(j),y(j),period,d2)
              if(d2 < r2) {
                # point j covers grid point
                covered = .true.
                goto 4242
              }
            }
          }
        }
        if(ixp1 <= n) {
          do j = ixp1,n {
            if(neigh(j)) {
              call dist2(xgrid,ygrid,x(j),y(j),period,d2)
              if(d2 < r2) {
                # point j covers grid point
                covered = .true.
                goto 4242
              }
            }
          }
        }
        if(.not. covered)
          kount = kount + 1
4242    continue
      }
    }
  }
}

# `kount' is the number of uncovered grid points

if(etainv > eps) {
  # usual calculation
  area = dble(kount) * dx * dy
  cifval = beta * exp(log(etainv) * area/(pi * r2))
} else {
  # etainv close to 0
  if(kount == 0)
    cifval = beta
  else
    cifval = zero
}

return
end
