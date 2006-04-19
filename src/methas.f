C Output from Public domain Ratfor, version 1.0
      subroutine methas(nmbr,aw,par,period,xprop,yprop,mprop,ntypes,ptyp
     *es,iseed, nrep,mrep,p,q,npmax,nverb,x,y,marks,aux,npts,fixall)
      implicit double precision(a-h,o-z)
      dimension par(1), ptypes(1), iseed(3)
      dimension x(1), y(1), marks(1), period(2)
      dimension xprop(1), yprop(1), mprop(1)
      integer aux(1)
      logical verb, marked, fixall
      eps = 2.22d-16
      one = 1.d0
      zero = 0.d0
      m1 = -1
      verb = .not.(nverb.eq.0)
      marked = ntypes .gt. 1
      irep = mrep
23000 if(irep .le. nrep)then
      if(verb .and. mod(irep,nverb).eq.0)then
      iprt = irep/nverb
      call intpr('irep/nverb=',-1,iprt,1)
      endif
      itype = 0
      call aru(1,zero,one,iseed,urp)
      if(urp.gt.p)then
      call aru(1,zero,one,iseed,urq)
      if(urq.gt.q)then
      u = xprop(irep)
      v = yprop(irep)
      if(marked)then
      mrk = mprop(irep)
      else
      mrk = -1
      endif
      call cif(nmbr,u,v,mrk,m1,x,y,marks,npts,ntypes, par,period,cifval,
     *aux)
      if(marked)then
      ffm = ptypes(mrk)
      else
      ffm = one
      endif
      a1 = q*aw*cifval/(ffm*(1-q)*(npts+1))
      a = min(one,a1)
      call aru(1,zero,one,iseed,bp)
      if(bp.lt.a)then
      npts = npts + 1
      itype = 1
      endif
      else
      if(npts.eq.0)then
      itype = 0
      else
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      if(marked)then
      mrki = marks(ix)
      ffm = ptypes(mrki)
      else
      mrki = -1
      ffm = one
      endif
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y, marks,npts,ntypes,par,perio
     *d,cifval,aux)
      a1 = (one-q)*npts*ffm/(q*aw*cifval)
      a = min(one,a1)
      call aru(1,zero,one,iseed,dp)
      if(dp.lt.a)then
      itype = 2
      endif
      endif
      endif
      else
      if(npts .eq. 0)then
      itype = 0
      else
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      u = xprop(irep)
      v = yprop(irep)
      if(marked)then
      mrki = marks(ix)
      if(fixall)then
      mrk = mrki
      else
      mrk = mprop(irep)
      endif
      ffm = ptypes(mrki)/ptypes(mrk)
      else
      mrki = -1
      mrk = -1
      ffm = one
      endif
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y,marks,npts, ntypes,par,perio
     *d,cvd,aux)
      call cif(nmbr,u,v,mrk,ix,x,y,marks,npts, ntypes,par,period,cvn,aux
     *)
      if(cvd .lt. eps)then
      if(cvn .lt. eps)then
      goto 23000
      else
      a = one
      endif
      else
      a1 = ffm*cvn/cvd
      a = min(one,a1)
      endif
      call aru(1,zero,one,iseed,sp)
      if(sp.lt.a)then
      itype = 3
      endif
      endif
      endif
      if(itype .gt. 0)then
      if(nmbr .eq. 8)then
      call updaux(itype,x,y,u,v,npts,ix,par,period,aux)
      endif
      if(itype.eq.1)then
      ix = npts
      x(ix) = u
      y(ix) = v
      if(marked)then
      marks(ix) = mrk
      endif
      if(npts .ge. npmax)then
      mrep = irep+1
      return
      endif
      else
      if(itype.eq.2)then
      call death(x,y,marks,marked,npts,ix)
      else
      x(ix) = u
      y(ix) = v
      if(marked)then
      marks(ix) = mrk
      endif
      endif
      endif
      endif
      irep = irep+1
      mrep = irep
      goto 23000
      endif
23001 continue
      return
      end
