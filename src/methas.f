      subroutine methas(nmbr,rw,par,period,tpar,ntypes,ptypes,iseed,
&     nrep,mrep, p,q,npmax,nverb,x,y,marks,aux,npts,fixall)
      implicit double precision(a-h,o-z)
      dimension par(1), tpar(1), ptypes(1), rw(4), iseed(3)
      dimension x(1), y(1), marks(1), period(2)
      integer aux(1)
      logical verb, marked, fixall
      eps = 2.22d-16
      one = 1.d0
      zero = 0.d0
      m1 = -1
      aw = (rw(2)-rw(1))*(rw(4)-rw(3))
      verb = .not.(nverb.eq.0)
      marked = ntypes .gt. 1
      irep = mrep
23000 if(.not.(irep .le. nrep))goto 23001
      if(.not.(verb .and. mod(irep,nverb).eq.0))goto 23002
      iprt = irep/nverb
      call intpr('irep/nverb=',-1,iprt,1)
23002 continue
      itype = 0
      call aru(1,zero,one,iseed,urp)
      if(.not.(urp.gt.p))goto 23004
      call aru(1,zero,one,iseed,urq)
      if(.not.(urq.gt.q))goto 23006
      call aru(1,rw(1),rw(2),iseed,u)
      call aru(1,rw(3),rw(4),iseed,v)
      if(.not.(marked))goto 23008
      call rmark(1,ntypes,ptypes,iseed,mrk)
      goto 23009
23008 continue
      mrk = -1
23009 continue
      call cif(nmbr,u,v,mrk,m1,x,y,marks,npts,ntypes, par,period,tpar,
&     cifval,aux)
      if(.not.(marked))goto 23010
      ffm = ptypes(mrk)
      goto 23011
23010 continue
      ffm = one
23011 continue
      a1 = q*aw*cifval/(ffm*(1-q)*(npts+1))
      a = min(one,a1)
      call aru(1,zero,one,iseed,bp)
      if(.not.(bp.lt.a))goto 23012
      npts = npts + 1
      if(.not.(npts .gt. npmax))goto 23014
      mrep = irep
      return
23014 continue
      itype = 1
23012 continue
      goto 23007
23006 continue
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      if(.not.(marked))goto 23016
      mrki = marks(ix)
      ffm = ptypes(mrki)
      goto 23017
23016 continue
      mrki = -1
      ffm = one
23017 continue
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y, marks,npts,ntypes,par,
&     period,tpar,cifval,aux)
      a1 = (one-q)*npts*ffm/(q*aw*cifval)
      a = min(one,a1)
      call aru(1,zero,one,iseed,dp)
      if(.not.(dp.lt.a))goto 23018
      itype = 2
23018 continue
23007 continue
      goto 23005
23004 continue
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      call aru(1,rw(1),rw(2),iseed,u)
      call aru(1,rw(3),rw(4),iseed,v)
      if(.not.(marked))goto 23020
      mrki = marks(ix)
      if(.not.(fixall))goto 23022
      mrk = mrki
      goto 23023
23022 continue
      call rmark(1,ntypes,ptypes,iseed,mrk)
23023 continue
      ffm = ptypes(mrki)/ptypes(mrk)
      goto 23021
23020 continue
      mrki = -1
      ffm = one
23021 continue
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y,marks,npts, ntypes,par,
&     period,tpar,cvd,aux)
      call cif(nmbr,u,v,mrk,ix,x,y,marks,npts, ntypes,par,period,tpar,
&     cvn,aux)
      if(.not.(cvd .lt. eps))goto 23024
      if(.not.(cvn .lt. eps))goto 23026
      goto 23000
23026 continue
      a = one
23027 continue
      goto 23025
23024 continue
      a1 = ffm*cvn/cvd
      a = min(one,a1)
23025 continue
      call aru(1,zero,one,iseed,sp)
      if(.not.(sp.lt.a))goto 23028
      itype = 3
23028 continue
23005 continue
      if(.not.(itype .gt. 0))goto 23030
      if(.not.(nmbr .eq. 8))goto 23032
      call updaux(itype,x,y,u,v,npts,ix,par,period,aux)
23032 continue
      if(.not.(itype.eq.1))goto 23034
      ix = npts
      x(ix) = u
      y(ix) = v
      if(.not.(marked))goto 23036
      marks(ix) = mrk
23036 continue
      goto 23035
23034 continue
      if(.not.(itype.eq.2))goto 23038
      call death(x,y,marks,marked,npts,ix)
      goto 23039
23038 continue
      x(ix) = u
      y(ix) = v
      if(.not.(marked))goto 23040
      marks(ix) = mrk
23040 continue
23039 continue
23035 continue
23030 continue
      irep = irep+1
      goto 23000
23001 continue
      return
      end
