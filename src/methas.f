C Output from Public domain Ratfor, version 1.0
      subroutine methas(nmbr,rw,par,period,tpar,nmarks,pmarks,iseed, nre
     *p,mrep,p,q,npmax,nverb,x,y,marks,aux,npts)
      implicit double precision(a-h,o-z)
      dimension par(1), tpar(1), pmarks(1), rw(4), iseed(3)
      dimension x(1), y(1), marks(1), period(2)
      integer aux(1)
      logical verb, marked
      eps = 2.22d-16
      one = 1.d0
      zero = 0.d0
      m1 = -1
      aw = (rw(2)-rw(1))*(rw(4)-rw(3))
      verb = .not.(nverb.eq.0)
      marked = nmarks .gt. 0
      if(mrep .eq. 0)then
      call aru(npts,rw(1),rw(2),iseed,x)
      call aru(npts,rw(3),rw(4),iseed,y)
      if(marked)then
      call rmark(npts,nmarks,pmarks,iseed,marks)
      endif
      call initaux(nmbr,par,period,x,y,npts,aux)
      mrep = 1
      endif
      irep = mrep
23004 if(irep .le. nrep)then
      if(verb .and. mod(irep,nverb).eq.0)then
      iprt = irep/nverb
      call intpr('irep/nverb=',-1,iprt,1)
      endif
      itype = 0
      call aru(1,zero,one,iseed,urp)
      if(urp.gt.p)then
      call aru(1,zero,one,iseed,urq)
      if(urq.gt.q)then
      call aru(1,rw(1),rw(2),iseed,u)
      call aru(1,rw(3),rw(4),iseed,v)
      if(marked)then
      call rmark(1,nmarks,pmarks,iseed,mrk)
      else
      mrk = -1
      endif
      call cif(nmbr,u,v,mrk,m1,x,y,marks,npts,nmarks, par,period,tpar,ci
     *fval,aux)
      if(marked)then
      ffm = pmarks(mrk)
      else
      ffm = one
      endif
      a1 = q*aw*cifval/(ffm*(1-q)*(npts+1))
      a = min(one,a1)
      call aru(1,zero,one,iseed,bp)
      if(bp.lt.a)then
      npts = npts + 1
      if(npts .gt. npmax)then
      mrep = irep
      return
      endif
      itype = 1
      endif
      else
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      if(marked)then
      mrki = marks(ix)
      ffm = pmarks(mrki)
      else
      mrki = -1
      ffm = one
      endif
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y, marks,npts,nmarks,par,perio
     *d,tpar,cifval,aux)
      a1 = (one-q)*npts*ffm/(q*aw*cifval)
      a = min(one,a1)
      call aru(1,zero,one,iseed,dp)
      if(dp.lt.a)then
      itype = 2
      endif
      endif
      else
      call aru(1,zero,one,iseed,xi)
      ix = 1 + int(npts*xi)
      call aru(1,rw(1),rw(2),iseed,u)
      call aru(1,rw(3),rw(4),iseed,v)
      if(marked)then
      call rmark(1,nmarks,pmarks,iseed,mrk)
      mrki = marks(ix)
      ffm = pmarks(mrki)/pmarks(mrk)
      else
      mrki = -1
      ffm = one
      endif
      call cif(nmbr,x(ix),y(ix),mrki,ix,x,y,marks,npts, nmarks,par,perio
     *d,tpar,cvd,aux)
      call cif(nmbr,u,v,mrk,ix,x,y,marks,npts, nmarks,par,period,tpar,cv
     *n,aux)
      if(cvd .lt. eps)then
      if(cvn .lt. eps)then
      goto 23004
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
      if(itype .gt. 0)then
      if(nmbr .eq. 8)then
      call updaux(itype,x,y,u,v,npts,ix,par,period,aux)
      endif
      if(itype.eq.1)then
      ix=npts
      x(ix) = u
      y(ix) = v
      if(marked)then
      marks(ix) = mrk
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
      goto 23004
      endif
23005 continue
      return
      end
